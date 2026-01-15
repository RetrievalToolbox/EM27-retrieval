function read_PROFFAST_spectrum(
        fname::String;
        cut_wavenumber::Union{Nothing, Tuple{<:Real, <:Real}}=nothing
    )

    # Create the output dictionary
    d = Dict()
    
    # Read sections, which are delimited by the '$' character
    open(fname, "r") do io

        # First section just contains descriptions
        s = String(readuntil(io, '$'))
        # Second section contains site info etc.
        s = String(readuntil(io, '$'))
        ss = split(s)
        d["sitename"] = String(ss[1])
        
        #=
            Construct date time from Date + fractional UT hours
        =#

        year = 2000 + parse(Int, ss[2][1:2])
        month = parse(Int, ss[2][3:4])
        day = parse(Int, ss[2][5:6])

        hour_frac = parse(Float64, ss[3])
        hour = Int(floor(hour_frac))

        minute_frac = (hour_frac - hour) * 60.0
        minute = Int(floor(minute_frac))

        second_frac = (minute_frac - minute) * 60
        second = Int(floor(second_frac))

        millisecond_frac = (second_frac - second) * 1000.0
        millisecond = Int(floor(millisecond_frac))
        
        d["datetime"] = DateTime(year, month, day, hour, minute, second, millisecond)
        # Julian date can be useful
        d["juliandate"] = datetime2julian(d["datetime"])
        #=
            Date portion end
        =#

        d["elevation_angle"] = parse(Float64, ss[4])
        d["viewing_azimuth"] = parse(Float64, ss[5])
        d["duration"] = parse(Float64, ss[6])
        d["latitude"] = parse(Float64, ss[7])
        d["longitude"] = parse(Float64, ss[8])
        d["altitude"] = parse(Float64, ss[9])
        d["spectrum_name"] = String(ss[10])

        # This third section contains "Filter, OPDmax [cm], semi FOV [rad]"
        s = String(readuntil(io, '$'))
        ss = split(s)
        d["filter"] = parse(Int, ss[1])
        d["OPD_max"] = parse(Float64, ss[2])
        d["semi_FOV"] = parse(Float64, ss[3])

        # This third section contains ILS (simple = 1, extended = 2)
        s = String(readuntil(io, '$'))
        ss = split(s)
        d["ILS"] = parse(Int, ss[1])

        # Now for modulation efficiency, separated by `,`, but there might be
        # spaces in there, in which case the split(s) function already splits the
        # two numbers..
        if ss[2][end] == ','
            d["modulation_efficiency"] = [parse(Float64, ss[2][1:end-1]), parse(Float64, ss[3])]
        else
            d["modulation_efficiency"] = [parse(Float64, x) for x in split(ss[2], ",")]
        end
        
        
        

        # This fourth section contains the start and end wavenumbers as well as the
        # spacing and total number of spectral samples
        s = String(readuntil(io, '$'))
        ss = split(s)
        d["start_wavenumber"] = parse(Float64, ss[1])
        d["end_wavenumber"] = parse(Float64, ss[2])
        d["delta_wavenumber"] = parse(Float64, ss[3])
        d["N_wavenumber"] = parse(Int, ss[4])

        # Generate the actual wavenumber array
        d["wavenumber"] = collect(LinRange(
                d["start_wavenumber"],
                d["end_wavenumber"],
                d["N_wavenumber"]
                )
            )
            
        # Reading the header is now done, and the rest of the file is mostly the spectrum,
        # with some empy characters in between
        s = String(readuntil(io, '$'))
        read(io, UInt8)
        read(io, UInt8)
        read(io, Float64)
        read(io, Float64)
        read(io, Int32)

        spec = zeros(d["N_wavenumber"])
        for k in 1:d["N_wavenumber"]
            spec[k] = read(io, Float32)
        end

        d["spectrum"] = spec
        
    end

    # User might want to cut down on the wavenumber range due to keep data volume small..
    if !isnothing(cut_wavenumber)
        cutout(d, cut_wavenumber, make_copy=false)
    end
    

    return d

end

function cutout(
        d_in::Dict,
        cut_wavenumber;
        make_copy=true
        )

    # Produce a deepcopy
    if make_copy
        d = deepcopy(d_in)
    else
        d = d_in
    end

    @assert cut_wavenumber[1] < cut_wavenumber[2]
    @assert cut_wavenumber[1] > d["wavenumber"][1]
    @assert cut_wavenumber[2] < d["wavenumber"][end]
    
    wn_select = searchsortedfirst.(Ref(d["wavenumber"]), cut_wavenumber)

    d["wavenumber"] = d["wavenumber"][wn_select[1]:wn_select[2]]
    d["start_wavenumber"] = d["wavenumber"][1]
    d["end_wavenumber"] = d["wavenumber"][end]
    d["N_wavenumber"] = length(d["wavenumber"])
    d["spectrum"] = d["spectrum"][wn_select[1]:wn_select[2]]

    if make_copy
        return d
    else
        return nothing
    end

end


function calculate_PROFFAST_ISRF(
        ME,
        PE,
        ν,
        Δν;
        N_ifg::Integer=80,
        apolin_calbias::Real=0.017,
        apoeff::Real=0.0,
        apophas::Real=0.0,
        N_izf::Integer=32,
        OPD_max::Real=1.8
        )

    # How many half-samples for ILS? (n_ilsgrid_ms * izf)
    N_ils_radius = 40 * N_izf

    modulat = zeros(2, N_ifg)
    
    a_val = zeros(2, N_ifg - 1)
    b_val = zeros(2, N_ifg - 1)
    
    Δν_ILS = Δν / N_izf

    # Output arrays
    ILS = zeros(N_ils_radius*2 + 1)
    ILS_wn = zeros(N_ils_radius*2 + 1)

    apolin = ME + apolin_calbias
    align_scale = ν / 7200.0
    OPD = zeros(N_ifg)

    for i in 1:N_ifg
    
        OPD_rel = Float64(i - 1) / Float64(N_ifg - 1)
        OPD[i] = OPD_rel * OPD_max + 1.0e-6
        
        # Self-apodization, NBM apodization
        term = (1.0 - OPD_rel * OPD_rel)
        modNBM = 0.152442 - 0.136176 * term + 0.983734 * term * term
        
        # Linear modulation loss
        modlin = (1.0 + (align_scale^1.5) * (apolin - 1.0) * OPD_rel)
        # Quadratic modulation loss
        modquad = (1.0 + align_scale * align_scale * apoeff * OPD_rel * OPD_rel)
        
        modulat[1,i] = modNBM * modlin * modquad
        modulat[2,i] = modulat[1,i] * tan(align_scale * apophas)
    
    end
    
    for i in 1:N_ifg - 1
    
        ΔOPD = 1.0 / (OPD[i+1] - OPD[i])
    
        for j in [1,2]
            a_val[j,i] = (modulat[j,i] * OPD[i+1] - modulat[j,i+1] * OPD[i]) * ΔOPD
            b_val[j,i] = (modulat[j,i+1] - modulat[j,i]) * ΔOPD
        end
    end
    
    for i in 0:N_ils_radius
    
        Δy = 2 * pi * (Δν_ILS * i + 1.0e-6)
        
        sum_g = 0.0
        sum_u = 0.0
        
        for j in 1:N_ifg - 1
            
            sinuy = sin(OPD[j+1] * Δy)
            sinly = sin(OPD[j] * Δy)
            cosuy = cos(OPD[j+1] * Δy)
            cosly = cos(OPD[j] * Δy)
            
            t1g = a_val[1,j] * (sinuy - sinly) / Δy
            t2g = b_val[1,j] * (
                cosuy + OPD[j+1] * Δy * sinuy
                - cosly - OPD[j] * Δy * sinly
                ) / (Δy * Δy)
    
            t1u = -a_val[2,j] * (cosuy - cosly) / Δy
            t2u = -b_val[2,j] * (
                sinly + OPD[j+1] * Δy * cosuy
                - sinuy - OPD[j] * Δy * cosly
                ) / (Δy * Δy)
    
            sum_g += t1g + t2g
            sum_u += t1u + t2u
            
        end
    
        ILS[1 + N_ils_radius - i] = sqrt(2 / pi) * (sum_g + sum_u)
        ILS[1 + N_ils_radius + i] = sqrt(2 / pi) * (sum_g - sum_u)
        
        ILS_wn[1 + N_ils_radius - i] = -Δν_ILS * i
        ILS_wn[1 + N_ils_radius + i] = Δν_ILS * i
    end;
    
    # Dampen ILS towards the rim
    for i in 1:N_ils_radius
        cosapo = cos(0.5 * pi * i * i / (N_ils_radius^2))
        ILS[1 + N_ils_radius - i] *= cosapo
        ILS[1 + N_ils_radius + i] *= cosapo
    end
    
    ILS ./= sum(ILS)
    
    return ILS_wn, ILS

end


function read_map_file(
    map_fname::String;
    convert_to_dry=true
    )

    
    # Open up the file and read where the actual atmosphere data is
    fid = open(map_fname, "r")
    
    row_header, row_unit = parse.(Ref(Int), split(readline(fid)))

    for l in 1:row_unit - 2
        global tmp = readline(fid)
    end
    unit_str = split(tmp, ",")

    close(fid)
    
    # Read in the contents of the map file (without units first)
    map_df = CSV.read(map_fname, DataFrame, delim=",", skipto=row_unit+1, header=row_header-1);

    # Sort by pressure (TOA at the top)
    sort!(map_df, "Pressure")
    
    # Turn the unit information into Unitful-parsable ones
    unitful_str = String[]
    for s in unit_str
        this_s = replace(s,  
            "_" => "/", 
            r"[0-9]+" => x -> "^"*x,
            "molecules" => "mol",
            "parts" => "1"
            )    
        push!(unitful_str, this_s)
    end

    # Now produce a Dictionary that contains unit-attached quantities
    out_dict = Dict(
        x => map_df[!, x] * uparse(unitful_str[i]) for (i,x) in enumerate(names(map_df))
            )

    if convert_to_dry

        # Calculate wet H2O VMR
        h2o_dry = (1 ./ map_df[!, "h2o"] .- 1).^(-1)
        map_df[!, "h2o"] .= h2o_dry

        # And then use that for all other gases
        for gas in ["hdo", "co2", "n2o", "co", "ch4", "hf", "o2"]
            map_df[!, gas] .= map_df[!, gas] .* (1 .+ h2o_dry)
        end
        
    end

    
    return out_dict
    
end