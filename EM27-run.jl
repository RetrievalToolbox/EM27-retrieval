using ArgParse
using CSV, DataFrames
using Dates
using Glob
using Interpolations
using LinearAlgebra
using MPI
using NetCDF
using Plots
using Printf
using Unitful
using Statistics 
using StatsBase
using YAML

using Pkg
Pkg.activate("../G3RT.jl")
using G3RT


using Profile
# Get rid of this
using Plots


# Load code to deal with PROFFAST pre-processing outputs
include("PROFFAST_tools.jl");
# Load code to deal with retrieval preparation
include("EM27-retrieval.jl");
# Load code that defines the forward model
include("forward_model.jl");
# Load code to process (retrieve) spectra
include("process_scene.jl");
# Load code to store results
include("results.jl")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--folder_spectra"
            help = "Path to folder with PROFFAST-generated spectra."
            arg_type = String
            required = true
        "--folder_map"
            help = "Path to folder with MAP atmospheres"
            arg_type = String
            required = true
        "--folder_pressure"
            help = "Path to folder with surface pressure data"
            arg_type = String
            required = true
        "--windows"
            help = "Which spectral windows to process (integers, separated by spaces)"
            arg_type = Int
            nargs = '+'
            required = true
        "--windows_config"
            help = "Path to YAML spectral window configuration"
            arg_type = String
            required = true
        "--output_basename"
            help = "Basename that will prepend the output CSV file."
            arg_type = String
            default = "output"
            required = false
    end

    return parse_args(s)
end

function main()

    #= 
        Set the float type used several components of the algorithm.
    
        So far, there seems to be no performance gain in using Float32 over
        Float64, of course the memory footprint is potentially smaller.
    =# 
    FT = Float64
    
    # Parse the command line arguments into a key => value dictionary
    parsed_args = parse_commandline()


    #=
        Count the occurences of "windows" passed by the user
    =#
    wmap = countmap(parsed_args["windows"])

    # Check if the user supplied duplicate windows
    if any(values(wmap) .> 1)

        dup_wins = [k for (k,v) in wmap if v > 1]
        error("Requested windows must be unique! Duplicates detected: $(dup_wins)")

    end

    @info "Reading in window configuration file $(parsed_args["windows_config"])"
    win_yml = YAML.load_file(parsed_args["windows_config"])

    # Grab the window keys from the configuration file
    # (to be used to check against the user window request)
    valid_windows = [convert(Int, x) for x in keys(win_yml)]
    sort!(valid_windows)

    # Check if the requested window is present in the config file
    for idx_want in parsed_args["windows"]
        if !(idx_want in valid_windows)
            error("Spectral window index $(idx_want) not found in " *
                "window configuration file: $(valid_windows)")
        end
    end

    #=
        Read in the measurements from the folder

        These are PROFFAST-preprocessed spectra, NOT the raw
        interferograms stored by OPUS.
    =#

    @info "Reading measurment files"
    measurements = read_measurements(
        parsed_args["folder_spectra"],
        win_yml,
        parsed_args["windows"]
    )

    #=
        Read in the measured surface pressure and
        creates a linear interpolation object so we can
        query the data at the times of measurement
    =#

    @info "Reading pressure files"
    pressure_df = read_pressure_files(
        parsed_args["folder_pressure"]
    )
    @info "Pressure data has $(nrow(pressure_df)) entries."

    # Create interpolation object, so we can calculate the surface pressure
    # at any given (Unix)time
    psurf_interp = linear_interpolation(
        pressure_df.UnixTime,
        pressure_df.BaroTHB40 * u"hPa"
    )

    #=
        Generate the spectral window objects
        from the list of user requested window
        indices.
    =#

    @info "Creating spectral window objects"
    swin_dict = create_spectral_windows(
            win_yml,
            parsed_args["windows"],
            buffer=15.0
        )

    #=
        Create gas objects for each window and read
        the appropriate spectroscopy data
    =#
    
    gas_dict = create_gases(
        win_yml,
        20,
        FT,
        parsed_args["windows"]
    )

    #=
        Produce ISRFs for each window (window index -> ISRF table)

        Note that this routine is copied straight from the PROFFAST
        source code, and generates ISRFs based on the central wavenumber
        and the supplied modulation efficiencies.
    
    =#
    isrf_dict = Dict()
    for idx in keys(measurements)
        @info "Calculating ISRF tables for $(win_yml[idx]["name"])"
        # Check if all modulation efficiencies are the same!
        un_me = unique((x -> x["modulation_efficiency"]).(values(measurements[idx])))

        if length(un_me) != 1
            error("Modulation efficiencies in measurements are not the same!")
        end
    
        fname = first(keys(measurements[idx]))
        s = measurements[idx][fname]
        
        # Assume that the first measurement here contains the modulation efficiency
        isrf = create_isrf_tables(s["wavenumber"], s["modulation_efficiency"]...)
        isrf_dict[idx] = isrf
    
    end

    # Interpolate MAP atmospheres to measurement times
    @info "Interpolating MAP files to measurment times"
    map_per_file = interpolate_map(
        parsed_args["folder_map"],
        measurements
        )

    # Create dispersion objects for each spectral window
    @info "Creating dispersion objects"
    dispersions = Dict(
        idx => create_dispersion(first(values(measurements[idx])), swin_dict[idx]) 
            for (idx,swin) in swin_dict
        )

    @info "Creating state vectors"
    sv_dict = Dict(idx =>
        create_statevector(
            Dict(dispersions[idx].spectral_window => dispersions[idx]),
            gas_dict[idx]
            )
        for idx in keys(swin_dict)
        )

    # Spectral window loop, each window is processed separately
    for swin_idx in sort!(collect(keys(swin_dict)))

        # This big dictonary will hold all the results
        all_results = Dict()

        if swin_idx in keys(all_results)
            swin_results = all_results[swin_idx]
        else
            all_results[swin_idx] = Dict()
            swin_results = all_results[swin_idx]
        end

        @info "Processing $(swin_dict[swin_idx].window_name)"
    
        #=
            Create buffers
            First, we have to estimate how big the buffer needs to be
        =#

        #= 
            How many total high-res spectral points do we need to do calculations?
            (this depends on the number of spectral windows, their respective width,
            and buffered by the width of the instrument line shape)
        =#
        
        N1 = Int(
            ceil(((x -> x.N_hires).(values(swin_dict)) |> sum) * 1.1)
            )
        
        # How many total spectral points do we have at instrument-level
        N2 = Int(ceil(first(values(measurements[swin_idx]))["N_wavenumber"] * 1.1))

        @info "Creating buffers"
        
        buf, oe_buf = create_buffers(
            FT,
            Dict(swin_dict[swin_idx] => dispersions[swin_idx]),
            gas_dict[swin_idx],
            sv_dict[swin_idx],
            N1,
            N2,
            "aux_data/hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc", # Solar model filename
            N_RT_lev=20 #length(first(values(map_per_file))["Pressure"])
        );

        # Set observer
        buf.scene.observer = G3RT.UplookingGroundObserver();

        #= 
        
            Loop over all measurements to perform the retrieval

        =#

        # (get a sorted list first, so we process in order of time)
        srt_fname = collect(keys(measurements[swin_idx]))
        sort!(srt_fname)
        
        for fname in srt_fname
            
            this_meas = measurements[swin_idx][fname]

            @info "Processing $(fname)"
            
            #=
                Insert location
            =#
            
            loc = G3RT.EarthLocation(
                this_meas["longitude"] * u"°",
                this_meas["latitude"] * u"°",
                this_meas["altitude"] * u"km"
                );
            
            buf.scene.location = loc;

            # Get the dispersion 
            d0 = this_meas["wavenumber"][1]
            d1 = this_meas["wavenumber"][2] - this_meas["wavenumber"][1]

            # Set the measurement dictonary "dispersion => spectrum"
            meas_dict = Dict(dispersions[swin_idx] => this_meas["spectrum"])
            # Set the noise standard deviation (white noise for FT)
            noise_std = this_meas["noise_std"]
            
            # Grab the atmosphere for this
            this_atm = map_per_file[fname]

            # Take the surface pressure from the interpolation object for this particular
            # time of measurement.
            this_psurf = psurf_interp(datetime2unix(this_meas["datetime"]))

            s = @time process_spectrum(
                meas_dict,
                this_meas["datetime"],
                90.0 - this_meas["elevation_angle"],
                this_psurf,
                [ # Dispersion coefficients initial
                    d0 - d1,
                    d1,
                    0.0
                ],
                Dict(swin_dict[swin_idx] => dispersions[swin_idx]),
                noise_std,
                this_atm,
                sv_dict[swin_idx],
                Dict(dispersions[swin_idx] => isrf_dict[swin_idx]),
                buf,
                oe_buf,
                max_iter=20,
                N_RT_lev=20
                );


            # If the H2O band has been processed, we should take the corrected
            # water vapor (or humidity) profile and insert it here..
            
            
            
            #=
                Store the results in a container
            =#
            
            store_results!(
                swin_results,
                # This can be a filename, for example
                fname,
                s,
                buf
                )
            
            chi2 = G3RT.calculate_chi2(s)
            xgases = G3RT.calculate_xgas(buf.scene.atmosphere)

            swin = swin_dict[swin_idx]
            p = plot(layout=(2, 1), size=(1000, 650))
                           
            y = G3RT.get_measured(s, swin);
            y2 = buf.rt_buf.radiance.I[buf.rt_buf.indices[swin]];
        
            ww = dispersions[swin_idx].ww
            #ww = 1e7 ./ dispersions[idx].ww
            
            plot!(ww, y, label="EM27/SUN", linewidth=2.0)
            plot!(ww, y2, label="Fit", linewidth=1.0, linestyle=:dash, linecolor="black")
        
            rms = sqrt(mean((100 .* (y2 .- y) ./ maximum(y)).^2))
            
            plot!(ww, 100 .* (y2 .- y) ./ maximum(y), label="Residual", subplot=2)
            title!(@sprintf("RMS: %.4f %%", rms), subplot=2, titlefont=(10))
            
            
            xlabel!("Wavenumber [cm\$^{-1}\$]", subplot=2, font=(10, "FreeSerif"))
            savefig("plots/fit_$(split(fname, "/")[end])_$(swin_idx).pdf")
            
            
        end # End measurement / file loop
       
        output_df = convert_to_df(all_results)

        # We store only scalar values into the table CSV, so we select them here
        # (meaning we skip AKs, pressure weights, etc.)
        scalar_select = [i for (i,n) in enumerate(names(output_df[swin_idx])) 
            if !(output_df[swin_idx][1,n] isa AbstractVector)]

        scalar_df = output_df[swin_idx][!, scalar_select]
        # Show the table to the user
        println(scalar_df)

    
        # Store scalar value of data frame as text file
        out_filename = "$(parsed_args["output_basename"])_$(swin_idx).csv"
        CSV.write(out_filename, scalar_df)

    end # End spectral index loop

    # But also store everything as a nice NetCDF file
    # (detailed results, including AKs and pressure weights etc.)

    

    
end # End main


main()