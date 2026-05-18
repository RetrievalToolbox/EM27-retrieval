function load_config()

    # Load the YAML configuration file
    win_yml = YAML.load_file("windows.yml")

    # Get a list of ALL spectral window indices, turn them into Integers, and sort
    # (useful to have a list of ALL windows)
    all_idx = keys(win_yml) |> collect .|> Int |> sort

    return (win_yml, all_idx)

end

function plot_measurements()

    # Plots two measurements from the two detectors
    sm = read_PROFFAST_spectrum("example_data/spectra/250722_131526SM.BIN");
    sn = read_PROFFAST_spectrum("example_data/spectra/250722_131526SN.BIN");

    plot(sm["wavenumber"], sm["spectrum"], label="SM")
    plot!(sn["wavenumber"], sn["spectrum"], label="SN")
    xlabel!("Wavenumber cm-1")
    ylabel!("Signal")

end


function first_filenames(measurements)

    return Dict(
        idx => collect(keys(measurements[idx])) |> sort |> first
        for idx in all_idx # Get them for ALL spectral windows
    )

end

function plot_more_measurements(measurements)

    for idx in sort(collect(keys(win_yml)))

        meas = measurements[idx][fname_first[idx]]

        p = plot()
        plot!(p, meas["wavenumber"], meas["spectrum"], label=nothing, size=(1000, 400))
        xlabel!("Wavenumber [cm-1]")

        title!("Spectral window #$(idx): $(win_yml[idx]["name"])")
        display(p)

    end

end

function create_gases(win_yml, idx_want)
    return create_gases(
        win_yml, # Windows YAML file that defines our spectral ranges
        20, # Number of pressure levels in model gas
        Float64, # Number type (leave this always as Float64)
        idx_want # List of spectral window identifiers
    )
end

function create_isrf(measurements)

    # Produce ISRFs for each window (window index -> ISRF table)
    isrf_dict = Dict()

    for idx in keys(measurements)

        # We only need the first measurement, they should all be the same anyway..
        s = first(values(measurements[idx]))

        isrf = create_isrf_tables(s["wavenumber"], s["modulation_efficiency"]...)
        isrf_dict[idx] = isrf

    end

    return isrf_dict

end

function plot_isrf(measurements, isrf_dict)

    plot_idx = 1 # Plot the ISRF for this spectral position
    nu_at_idx = first(values(measurements[idx_want]))["wavenumber"][plot_idx] # Get the wavenumber

    p = plot()

    plot!(
        p,
        isrf_dict[idx_want].ww_delta[:, plot_idx],
        isrf_dict[idx_want].relative_response[:, plot_idx],
        lw=2.0,
        label="ν = $(nu_at_idx) cm-1"
    )

    xlabel!("Δν [cm\$^{-1}\$]\n")
    ylabel!("\nISRF [cm]")

end


function create_dispersions(measurements, swin_list)

    dispersions = Dict(
        idx => create_dispersion(
                first(values(measurements[idx])),
                swin
        ) for (idx, swin) in swin_list
    )

    for swin in values(swin_list)
        println("$(swin) has $(swin.N_hires) spectral points")
    end

    return dispersions

end

function plot_model_atmosphere(measurements)

    p = plot(size=(450, 550))
    for fname in keys(measurements[1])

        map_file = map_per_file[fname]

        plot!(
            p,
            map_file["h2o"],
            map_file["Pressure"],
            label=nothing,
            );
    end
    yflip!()

end


function plot_pressure_data(measurements)

    _plot_dates = DateTime[]
    _plot_pressure = Float64[]
    for (_, meas) in measurements[idx_want]
        push!(_plot_dates, meas["datetime"])
        push!(_plot_pressure, meas["pressure"] |> ustrip)
    end

    scatter(_plot_dates, _plot_pressure, label=nothing,
            markershape=:circle, size=(1000, 400))
    ylabel!("Surface pressure")
    xlabel!("Time of measurement UTC")

end


function create_state_vector(dispersions, swin_list)

    # A dictionary that maps spectral window => dictonary

    return Dict(idx =>
        create_statevector(
            Dict(dispersions[idx].spectral_window => dispersions[idx]),
            gas_dict[idx], # A vector of gas objects
            1, # Order of solar continuum polynomial
            1, # Order of dispersion polynomial
        )
        for idx in keys(swin_list)
    )

end


function create_buffer(measurements, swin_list, dispersions, gas_dict, sv_dict)

    N1 = swin_list[idx_want].N_hires * 1.1 |> ceil |> Int
    N2 = measurements[idx_want][fname_first[idx_want]]["N_wavenumber"] * 1.1 |> ceil |> Int

    # Create the `buffer` object!
    buf = create_buffers(
        Float64, # Number type to be used for the vectors and arrays
        Dict(d.spectral_window => d for (k,d) in dispersions), # Dictionary [spectral window => dispersion]
        gas_dict[idx_want], # List of gases
        sv_dict[idx_want], # State vector(s)
        N1,
        N2;
        N_RT_lev=20 # 20 levels for our RT calculations!
    )

    return buf

end


function ingest_location!(buf, measurements)

    loc = RE.EarthLocation(
        measurements[idx_want][fname_first[idx_want]]["longitude"],
        measurements[idx_want][fname_first[idx_want]]["latitude"],
        measurements[idx_want][fname_first[idx_want]]["altitude"],
        u"km"
    )

    buf.scene.location = loc
    buf.scene.observer = RE.UplookingGroundObserver()

    # Show the new location object with proper values
    buf.scene.location

end


function process_spectrum(file_idx, all_fnames, measurements, buf;
    max_iter=20
    )

    fname = all_fnames[file_idx] # File name corresponding to `file_idx`
    meas = measurements[idx_want][fname] # Measurement data dictionary corresponding to `fname`

    swin = buf.spectral_window[1]
    disp = buf.rt_buf.dispersion[swin]

    # Dispersion coefficients derived from the measurement
    d0 = meas["start_wavenumber"]
    d1 = meas["delta_wavenumber"]

    s = process_spectrum(
        Dict(disp => meas["spectrum"]), # Measured spectrum as dictonary: dispersion => spectrum
        meas["datetime"], # Datetime of measurement
        90.0 - meas["elevation_angle"], # Solar zenith angle
        meas["pressure"], # Measured surface pressure (with physical unit!)
        [ # Dispersion coefficients initial
            d0 - d1,
            d1,
            0.0
        ],
        Dict(disp.spectral_window => disp), # Spectral window => dispersion
        meas["noise_std"], # Estimated noise, expressed as the standard deviation of a Gaussian in signal units
        map_per_file[fname], # Model atmosphere prior for this measurement
        sv_dict[idx_want], # The state vector object
        Dict(dispersions[i] => isrf_dict[i] for i in keys(swin_list)), # Dispersion => instrument spectral response function
        buf; # The buffer object
        max_iter=max_iter # Number of iterations before we give up
    )

    return s

end

function plot_fit_result(s, buf)

    swin = buf.spectral_window[1]

    p = plot(layout=(2, 1), size=(1000, 650))
    snr_label = @sprintf "SNR: %.0f" ((RE.get_measured(s) ./ RE.get_noise(s)) |> mean)


    y = RE.get_measured(s, swin);
    y2 = RE.get_modeled(s, swin);
    eps = RE.get_noise(s, swin);

    ww = RE.get_wavenumber(s, swin);

    plot!(ww, y, label="EM27/SUN", linewidth=2.0) #, ylims=(0, 9))
    plot!(ww, y2, label="Fit", linewidth=1.0, linestyle=:dash, linecolor="black")
    plot!(ww, y2 .- y, label="Residual", linewidth=1.0, linecolor="black")

    title!(snr_label, subplot=1)
    ylabel!("Signal", subplot=1)
    xlabel!("Wavenumber [cm\$^{-1}\$]", subplot=2, font=(10, "FreeSerif"))

    rrms = sqrt(mean((100.0 * (y2 .- y) ./ y).^2))

    plot!(ww, (y2 .- y) ./ eps, label="Residual", subplot=2)
    title!(@sprintf("RRMS: %.4f %%", rrms), subplot=2, titlefont=(10))
    ylabel!("Fractions of noise", subplot=2)



    xlabel!("Wavenumber [cm\$^{-1}\$]", subplot=2, font=(10, "FreeSerif"))
    display(p)

end


function plot_latest_iteration(solver, buf)

    p = plot(layout=(2, length(swin_list)), size=(1000, 650))

    for (i,swin) in enumerate(values(swin_list))

        y = RE.get_measured(solver, swin);
        y2 = RE.get_modeled(solver, swin);
        eps = RE.get_noise(solver, swin);

        ww = RE.get_wavenumber(solver, swin);

        plot!(ww, y, label="EM27/SUN", linewidth=2.0) #, ylims=(0, 9))
        plot!(ww, y2, label="Fit", linewidth=1.0, linestyle=:dash, linecolor="black")

        title!("Iteration #$(RE.get_iteration_count(solver))", subplot=1, titlefont=(10))
        ylabel!("Signal", subplot=1)
        xlabel!("Wavenumber [cm\$^{-1}\$]", subplot=2, font=(10, "FreeSerif"))

        rrms = sqrt(mean((100.0 * (y2 .- y) ./ y).^2))

        plot!(ww, (y2 .- y) ./ eps, label="Residual", subplot=2)
        title!(@sprintf("RRMS: %.4f %%", rrms), subplot=2, titlefont=(10))
        ylabel!("Fractions of noise", subplot=2)

    end

    xlabel!("Wavenumber [cm\$^{-1}\$]", subplot=2, font=(10, "FreeSerif"))
    display(p)

    RE.print_state_vector_update(solver.state_vector);

end


function plot_timeseries(result_dict)

    # Load the results of the PROFFAST retrievals
    df_ref = CSV.File("./example_data/comb_invparms_GSFC_SN250_250722-250722.csv") |> DataFrame;
    # Subselect to the measurements we also have:
    df_specnames = df_ref.var" spectrum" .|> lstrip
    spec_names = (x -> split(x, "/")[end] |> String).(all_fnames)
    spec_names = (x -> replace(x, "SM" => "SN")).(spec_names) # Have to replace SM -> SN for this
    keep_idx = searchsortedfirst.(Ref(df_specnames), spec_names);

    df_ref = df_ref[keep_idx, :];
    # Make DateTimes from julian date
    df_ref[!, "Date"] .= julian2datetime.(df_ref.var" JulianDate");

    if "CO2" in keys(result_dict)
        p = plot()
        scatter!(p, result_dict["Date"], result_dict["CO2"] * 1e6, markershape=:circle, label="GSFC", size=(1000, 300))
        plot!(p, df_ref.Date, df_ref.var" XCO2_STR", markershape=:star, label="PROFFAST XCO2_STR")
        scatter!(p, df_ref.Date, df_ref.var" XCO2", markershape=:square, label="PROFFAST XCO2")
        #plot!(twinx(), result_dict["Date"], 1 ./ cosd.(result_dict["SZA"]), yaxis = "Airmass",
        #    color=:black, label=nothing, linewidth=2)
        title!("XCO2")
        display(p)
    end

    if "H2O" in keys(result_dict)
        p = plot()
        scatter!(p, result_dict["Date"], result_dict["H2O"] * 1e6, markershape=:circle, label="GSFC", size=(900, 300))
        scatter!(p, df_ref.Date, df_ref.var" XH2O", markershape=:square, label="PROFFAST")
        #plot!(twinx(), result_dict["Date"], 1 ./ cosd.(result_dict["SZA"]), yaxis = "Airmass",
        #    color=:black, label="Airmass", linewidth=2)
        title!("XH2O")
        display(p)
    end

    if "CH4" in keys(result_dict)
        p = plot()
        scatter!(p, result_dict["Date"], result_dict["CH4"] * 1e9, markershape=:circle, label="GSFC", size=(1000, 300))
        scatter!(p, df_ref.Date, df_ref.var" XCH4" * 1e3, markershape=:square, label="PROFFAST")
        scatter!(p, df_ref.Date, df_ref.var" XCH4_S5P" * 1e3, markershape=:square, label="PROFFAST")
        #plot!(twinx(), result_dict["Date"], 1 ./ cosd.(result_dict["SZA"]), yaxis = "Airmass",
        #    color=:black, label="Airmass", linewidth=2)
        title!("XCH4")
    end

    if "CO" in keys(result_dict)
        p = plot()
        scatter!(p, result_dict["Date"], result_dict["CO"] * 1e9, markershape=:circle, label="GSFC", size=(1000, 300))
        scatter!(p, df_ref.Date, df_ref.var" XCO" * 1e3, markershape=:square, label="PROFFAST")
        #plot!(twinx(), result_dict["Date"], 1 ./ cosd.(result_dict["SZA"]), yaxis = "Airmass",
        #    color=:black, label="Airmass", linewidth=2)
        title!("XCO")
    end

    if "O2" in keys(result_dict)
        plot(p, result_dict["Date"], result_dict["O2"], markershape=:circle, label="GSFC", size=(1000, 300))
        #plot!(twinx(), result_dict["Date"], 1 ./ cosd.(result_dict["SZA"]), yaxis = "Airmass",
        #    color=:green, label=nothing, linewidth=2)
        title!("XO2")
    end

end