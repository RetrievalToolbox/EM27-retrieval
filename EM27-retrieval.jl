"""
    Creates retrieval windows from the YAML-parsed dict.
"""
function create_spectral_windows(
        d::Dict,
        idx_want::Union{T, Vector{T}};
        buffer=15.0
    ) where {T <: Integer}

    swin_dict = Dict{Integer, RE.SpectralWindow}()

    # Make this a 1-element list if user supplied a number
    if idx_want isa Integer
        idx_want = [idx_want]
    end
    
    # Loops through all of the integer-based keys
    for idx in keys(d)
        
        # Only process windows the user wants us to process
        if !(idx in idx_want)        
            # @info "Skipping [$(idx)] / $(d[idx]["name"])"
            continue
        end

        @info "Adding $(d[idx]["name"]) to the list of retrieval windows."

        wn_start = convert(Float64, d[idx]["wavenumber_start"])
        wn_end = convert(Float64, d[idx]["wavenumber_end"])
        wn_ref = convert(Float64, d[idx]["wavenumber_reference"])

        # Create the high-resolution model grid
        ww = collect(wn_start - buffer:0.01:wn_end + buffer)
        
        swin_dict[idx] = RE.SpectralWindow(
            d[idx]["name"],
            wn_start,
            wn_end,
            ww,
            u"cm^-1",
            wn_ref
        )

    end

    return swin_dict

end


"""
    Creates dispersion objects, given a PROFFAST-type spectrum dict
"""
function create_dispersion(
        spec_dict::Dict,
        swin::RE.SpectralWindow
    )

    # Polynomial dispersion is straightforward due to the 
    # linear nature of the FTIR

    nu0 = spec_dict["wavenumber"][1] * u"cm^-1"
    delta_nu = (spec_dict["wavenumber"][2] - spec_dict["wavenumber"][1]) * u"cm^-1"
    
    return RE.SimplePolynomialDispersion(
        [nu0, delta_nu, 0.0u"cm^-1", 0.0u"cm^-1"],
        1:spec_dict["N_wavenumber"],
        swin
    )
    
end

function read_measurements(
        folder_spectra::String,
        d::Dict,
        idx_want::Union{T, Vector{T}};
        wn_buffer=25.0
) where {T <: Integer}

    if idx_want isa Integer
        idx_want = [idx_want]
    end

    # Find all spectra files
    flist_sm = glob(joinpath(folder_spectra, "*SM.BIN"))
    sort!(flist_sm)
    flist_sn = glob(joinpath(folder_spectra, "*SN.BIN"))
    sort!(flist_sn)

    # Read them all in (the whole contents)
    fdict = Dict()
    fdict["SM"] = Dict(
        fname => read_PROFFAST_spectrum(fname) for fname in flist_sm
        )
    
    fdict["SN"] = Dict(
        fname => read_PROFFAST_spectrum(fname) for fname in flist_sn
        )

    # Now we read the window dictionary `d` and assign a `filename` -> `dict`
    # to each.

    spectra_dict = Dict{Integer, Dict}()

    for idx in keys(d)
        
        # Only process windows the user wants us to process
        if !(idx in idx_want)        
            @info "Skipping [$(idx)] / $(d[idx]["name"])"
            continue
        else
            @info "Adding [$(idx)] / $(d[idx]["name"])"
        end

        cut_wavenumber = [
            d[idx]["wavenumber_start"] - wn_buffer,
            d[idx]["wavenumber_end"] + wn_buffer
            ]

        this_fdict = fdict[d[idx]["detector"]]
        
        subdict = Dict{String, Dict}()
        for fname in sort!(collect(keys(this_fdict)))

            # Before cutting, however, we grab the noise level from
            # the full spectral range
            
            noise = estimate_noise(this_fdict[fname])
                        
            subdict[fname] = cutout(
                this_fdict[fname],
                cut_wavenumber
            )
            subdict[fname]["noise_std"] = noise

        end

        spectra_dict[idx] = subdict
        
    end
    
    return spectra_dict
    
end

function estimate_noise(d::Dict)

    # First, we check the upper wavenumber bound
    if d["end_wavenumber"] > 10000
        noise_idx = searchsortedfirst.(Ref(d["wavenumber"]), 
            [12500, min(14000, d["end_wavenumber"])]
        )
    else
        noise_idx = searchsortedfirst.(Ref(d["wavenumber"]), 
            [max(3800, d["start_wavenumber"]), 3950]
        )
    end

    this_noise = std(@views d["spectrum"][:][noise_idx])

    return this_noise

end

"""
    Interpolates MAP files to the actual times of measurement
"""
function interpolate_map(
        folder_map::String,
        m::Dict
        )
        
    # List all MAP files
    map_flist = sort!(glob(joinpath(folder_map, "*.map")))
    
    # Extract time strings from filenames
    tstrings = (x -> x.match).(match.(Ref(r"[0-9]{10}"), map_flist))
    
    # Parse them into DateTimes
    datetimes = Dates.DateTime.(tstrings, Ref(dateformat"yyyymmddHH"))

    # We will now create a dictionary with keys that are equal to the
    # measurement times

    map_dict = Dict{String, Dict{String, Vector}}()
    
    
    for idx in keys(m)
        for fname in keys(m[idx])

            # Check if this exact file has already been processed..
            if fname in keys(map_dict)
                continue
            end

            dt = m[idx][fname]["datetime"]
            
            if length(map_flist) == 1
                idx_left = 1
                idx_right = 1
                fac = 0.5
            else
                idx_right = searchsortedfirst(datetimes, dt)
                idx_left = idx_right - 1
                fac = (dt - datetimes[idx_left]) / (datetimes[idx_right] - datetimes[idx_left])
            end
            
            

            map_left = read_map_file(map_flist[idx_left])
            map_right = read_map_file(map_flist[idx_right])

            this_map = map_left
            for k in keys(this_map)
                this_map[k] .*= fac
                this_map[k] .+= map_right[k] .* (1 - fac)
            end

            map_dict[fname] = this_map
            
        end
    end

    return map_dict
    
end

function copy_map_to_atmosphere!(
    scene::RE.EarthScene,
    map::Dict
    )

    atm = scene.atmosphere
    # `map` is a dictionary with Unitful arrays, so we can do unit conversions here!

    # MET pressure grid
    @views atm.met_pressure_levels[:] = ustrip.(Ref(atm.met_pressure_unit), copy(map["Pressure"]))
    @views atm.met_pressure_layers[:] = RE.levels_to_layers(atm.met_pressure_levels)

    # MET temperature grid
    @views atm.temperature_levels[:] = ustrip.(Ref(atm.temperature_unit), copy(map["Temp"]))
    @views atm.temperature_layers[:] = RE.levels_to_layers(atm.temperature_levels)

    # Calculate SH from H2O
    sh = RE.H2O_VMR_to_specific_humidity.(map["h2o"])
    @views atm.specific_humidity_levels[:] = ustrip.(Ref(atm.specific_humidity_unit), sh)
    @views atm.specific_humidity_layers[:] = RE.levels_to_layers(atm.specific_humidity_levels)

    # Compute altitude and gravity
    RE.calculate_altitude_and_gravity!(scene)

    # .. and mid-layer values
    atm.gravity_layers[:] = RE.levels_to_layers(atm.gravity_levels)
    atm.altitude_layers[:] = RE.levels_to_layers(atm.altitude_levels)

end

function create_statevector(
        dispersions::Dict{<:RE.AbstractSpectralWindow, <:RE.AbstractDispersion},
        gas_list::Vector{RE.GasAbsorber},
        solar_order::Integer,
        dispersion_order::Integer
)

    solar_array = RE.SolarScalerPolynomialSVE[]
    
    for swin in keys(dispersions)
        for p in 0:solar_order

            if p == 0
                fg = 1.0
            else
                fg = 0.0
            end
            my_solar_scale_sve = RE.SolarScalerPolynomialSVE(
                swin,
                p,
                u"cm^-1",
                fg,
                fg,
                10.0
            )

            push!(solar_array, my_solar_scale_sve)
            
        end
    end

    
    disp_array = RE.DispersionPolynomialSVE[]
    
    for disp in values(dispersions) # Loop over dispersion objects

        for p in 0:dispersion_order # Loop over polynomial order
        
            push!(disp_array, 
                RE.DispersionPolynomialSVE(
                    disp,
                    p, # Polynomial order
                    u"cm^-1", # Unit
                    0.0, # First guess
                    0.0, # Prior value
                    0.0, # Prior covariance
                )
            )
        
        end
    
    end

    zlo_array = RE.ZeroLevelOffsetPolynomialSVE[]
    for swin in keys(dispersions)
        
        zlo_sve0 = RE.ZeroLevelOffsetPolynomialSVE(
            swin, # Retrieval window name
            0,
            u"cm^-1", # wavelength units
            Unitful.NoUnits, # radiance units
            0.0, # First guess
            0.0, # Prior value
            (10.0)^2  # Prior cov
        );

        push!(zlo_array, zlo_sve0)

        zlo_sve1 = RE.ZeroLevelOffsetPolynomialSVE(
            swin, # Retrieval window name
            1,
            u"cm^-1", # wavelength units
            Unitful.NoUnits, # radiance units
            0.0, # First guess
            0.0, # Prior value
            (10.0)^2  # Prior cov
        );

        push!(zlo_array, zlo_sve1)
    end
    
    gas_scale_array = RE.GasLevelScalingFactorSVE[]
    for gas in gas_list
        
        this_sve = RE.GasLevelScalingFactorSVE(
            1,
            length(gas.vmr_levels),
            gas,
            Unitful.NoUnits,
            1.0,
            1.0,
            1.0
            )

        push!(gas_scale_array, this_sve)
        
    end

    my_temp_sve = RE.TemperatureOffsetSVE(
        u"K",
        0.0,
        0.0,
        (5.0)^2
        )

    my_psurf_sve = RE.SurfacePressureSVE(
        u"hPa",
        0.0,
        0.0,
        (25.0)^2
        )

    return RE.RetrievalStateVector([
            solar_array...,
            #my_temp_sve,
            #my_psurf_sve,
            disp_array...,
            gas_scale_array...,
            #zlo_array...,
            ])

end

function create_gases(
    d::Dict,
    N_levels::Integer,
    my_type::Type{T},
    idx_want::Union{U, Vector{U}}
) where {T, U <: Integer}

    gas_dict = Dict{Integer, Vector{RE.GasAbsorber}}()

    if idx_want isa Integer
        idx_want = [idx_want]
    end

    for idx in keys(d)

        # Only process windows the user wants us to process
        if !(idx in idx_want)        
            # @info "Skipping $(d[idx]["name"])"
            continue
        end
        
        gas_sublist = RE.GasAbsorber[]
        
        for spec in keys(d[idx]["spectroscopy"])

            @info "Loading $(d[idx]["spectroscopy"][spec]["ABSCO"])"
        
            if haskey(d[idx]["spectroscopy"][spec], "scale_factor")
                scale_factor = convert(
                    Float64,
                    d[idx]["spectroscopy"][spec]["scale_factor"]
                )
            else
                scale_factor = 1.0
            end

            # Read in the ABSCO
            absco = RE.load_ABSCOAER_spectroscopy(
                d[idx]["spectroscopy"][spec]["ABSCO"],
                spectral_unit=:Wavenumber,
                scale_factor=scale_factor
            )

            # parse the units
            this_unit_str = d[idx]["spectroscopy"][spec]["unit"]
            if this_unit_str == "1"
                this_unit = Unitful.NoUnits
            else
                this_unit = Unitful.uparse(this_unit_str)
                if !(this_unit isa Unitful.DimensionlessUnits)
                    error("VMR unit must dimensionless: $(this_unit)")
                end
            end
            
            # Create the gas object with correct units
            gas = RE.GasAbsorber(
                d[idx]["spectroscopy"][spec]["name"],
                absco,
                zeros(my_type, N_levels),
                this_unit
            )

            # Push this gas object into the sublist
            # (list of gas objects for this spectral window)
            push!(gas_sublist, gas)

        end

        # Add this gas sublist into the dictonary
        gas_dict[idx] = gas_sublist

    end

    return gas_dict

end


function create_buffers(
        my_type, 
        dispersions, 
        gas_list,
        sv, 
        N1, # How many spectral points needed including convolution?
        N2; # How many points in total needed for instrument-level radiances?
        N_RT_lev=20
    )

    # Read in the solar model, for each spectral window
    #solar_model_dict = Dict(swin => RE.TSISSolarModel(
    #        solar_model_fname,
    #        spectral_unit=:Wavenumber
    #    ) for swin in keys(dispersions)
    #);

    solar_fname = "./aux_data/solar_merged_20240731_600_33300_000.out"
    solar_df = CSV.File(
        solar_fname,
        delim=' ',
        ignorerepeated=true,
        skipto=4
    ) |> DataFrame

    # Only select 2000cm-1 through 9000cm-1 (like PROFFAST does)
    solar_select = searchsortedfirst.(Ref(solar_df[:,1]), [2000., 9000.])
        
    solar_model_dict = Dict(
        swin => RE.ListSolarModel(
            solar_df[solar_select[1]:solar_select[2], 1],
            u"cm^-1",
            solar_df[solar_select[1]:solar_select[2], 2],
            Unitful.NoUnits
        ) for swin in keys(dispersions)
    );

    # In order for this to work in arbitrary solar irradiance units, we
    # must set the units manually to 1
    #for solar_model in values(solar_model_dict)
    #    solar_model.irradiance_unit = Unitful.NoUnits
    #    max_irrad = maximum(solar_model.irradiance)
    #    @views solar_model.irradiance[:] ./= max_irrad
    #end
    
    # Will contain outputs of ISRF application
    inst_buf = RE.InstrumentBuffer(
        zeros(my_type, N1),
        zeros(my_type, N1),
        zeros(my_type, N2),
    );
    

    # Will contain outputs of radiance, jacobians and
    # dispersion indices
    
    rt_buf = RE.ScalarRTBuffer(
        dispersions, # Already a SpectralWindows -> Dispersion dictionary
        RE.ScalarRadiance(my_type, N2), # Hold the radiance
        Dict(sve => RE.ScalarRadiance(my_type, N2) for sve in sv.state_vector_elements),
        Dict(swin => zeros(Int, 0) for swin in keys(dispersions)), # Hold the detector indices
        Unitful.NoUnits, # Radiance units for the forward model, we use a DOAS-type approach for EM-27 where we normalize
    );

    # Just grab a list of the solar models
    swin_list = collect(keys(dispersions))
    
    # Finally, create the big containter which links them all

    #=
    buf = RE.EarthAtmosphereBuffer(
        sv, # The state vector
        swin_list, # The spectral window (or a list of them)
        [(:NoSurface, 0) for x in swin_list], # Surfaces (none needed for uplooking)
        [gas_list...,], # All atmospheric elements
        solar_model_dict, # Solar model dictionary (spectral window -> solar model)
        [:XRTM for swin in swin_list], #
        rt_buf,
        inst_buf,
        N_RT_lev, # The number of retrieval or RT pressure levels
        50, # The number of meteorological pressure levels, as given by the atmospheric inputs
        my_type # The chosen Float data type (e.g. Float16, Float32, Float64)
    );
    =#
    buf = RE.EarthAtmosphereBuffer(
        sv, # The state vector
        swin_list, # The spectral window (or a list of them)
        [(:NoSurface, 0) for x in swin_list], # Surfaces (none needed for uplooking)
        [gas_list...,], # All atmospheric elements
        solar_model_dict, # Solar model dictionary (spectral window -> solar model)
        [:BeerLambert for swin in swin_list], #
        RE.ScalarRadiance,
        rt_buf,
        inst_buf,
        N_RT_lev, # The number of retrieval or RT pressure levels
        50, # The number of meteorological pressure levels, as given by the atmospheric inputs
        my_type # The chosen Float data type (e.g. Float16, Float32, Float64)
    );

    return buf

end


function create_isrf_tables(
        nu_grid,
        ME::Real, # Modulation efficiency
        PE::Real # Phase efficiency
    )

    N_nu = length(nu_grid)
    N_ifg = 80
    N_izf = 32
    N_delta_nu = N_ifg * N_izf + 1
    
    isrf_delta_nu = zeros(Float32, N_delta_nu, N_nu)
    isrf_rr = zeros(Float32, N_delta_nu, N_nu)

    for (i,nu) in enumerate(nu_grid)

        if i == N_nu
            d_nu = nu_grid[end] - nu_grid[end-1]
        else
            d_nu = nu_grid[i+1] - nu_grid[i]
        end
        
        v1, v2 = calculate_PROFFAST_ISRF(
            ME,
            PE,
            nu,
            d_nu,
            N_ifg=N_ifg,
            N_izf=N_izf
            )

        isrf_delta_nu[:,i] = v1
        isrf_rr[:,i] = v2
        
        
    end

    return RE.TableISRF(
        isrf_delta_nu,
        u"cm^-1",
        isrf_rr
        )
    
    
end



"""
    Reads pressure measurements from a csv file and adds the 
    relevant surface pressure data into the `measurements` dictionary. Uses the
    closest measurement time for now.
"""
function interpolate_pressure_files_B33!(measurements::Dict, fname::String)

    # Produce a concatenated DataFrame containing ALL pressure data
    df = CSV.File(fname) |> DataFrame
    
    
    # Fill in the correct date/times from TIME_UTC
    df[!, "FullDate"] .= (x -> Dates.DateTime(x, Dates.dateformat"yyyy-mm-dd HH:MM:SS")).(df.TIME_UTC)

    # Now loop over all spectral windows:
    for idx in keys(measurements)
        for fname in keys(measurements[idx])

            # Grab this measurement belonging to spectral window idx and filename fname
            meas = measurements[idx][fname]
            
            # Find the time step in the pressure data which is closest to the 
            # EM27 measurement time.
            p_idx = searchsortedfirst(df.FullDate, meas["datetime"])

            # Add the pressure data to the measruement dictionary!
            meas["pressure"] = df.BP_mbar[p_idx] * u"hPa"

        end
    end

end

"""
    Reads pressure measurements from a folder with `*.dat` files and adds the 
    relevant surface pressure data into the `measurements` dictionary. Uses the
    closest measurement time for now.
"""
function interpolate_pressure_files!(measurements::Dict, folder::String)

    fnames = glob(joinpath(folder, "*.dat")) |> sort

    # Produce a concatenated DataFrame containing ALL pressure data
    df = vcat((CSV.File.(fnames) .|> DataFrame)...)
    
    # Create new column with DateTimes
    df[!, :FullDate] .= Dates.DateTime(1999, 1, 1, 0,0,0)
    
    # Fill in the correct date/times from UTCdate and UTCtime
    for row in eachrow(df)
        this_date_str = row.UTCdate_____
        this_time_str = row.UTCtime___ |> string
    
        this_fulldatetime = Dates.DateTime("$(this_date_str) $(this_time_str)", dateformat"dd.mm.yyyy HH:MM:SS")
        
        row.FullDate = this_fulldatetime
    end

    # Now loop over all spectral windows:
    for idx in keys(measurements)
        for fname in keys(measurements[idx])

            # Grab this measurement belonging to spectral window idx and filename fname
            meas = measurements[idx][fname]
            
            # Find the time step in the pressure data which is closest to the 
            # EM27 measurement time.
            p_idx = searchsortedfirst(df.FullDate, meas["datetime"])

            # Add the pressure data to the measruement dictionary!
            meas["pressure"] = df.BaroTHB40[p_idx] * u"hPa"

        end
    end

end
