"""
    Grabs results from the various containers (buffer, solver) and stores them
    in an easy-to-read dictionary.
"""
function store_results!(
    in_dict::Dict,
    # This can be a filename, for example
    store_key,
    s::G3RT.AbstractSolver,
    buf::G3RT.EarthAtmosphereBuffer
    )

    d = Dict()
    in_dict[store_key] = d

    # Grab the only spectral window we have
    swin = first(keys(buf.rt))

    # Store the surface pressure used
    d["psurf"] = buf.scene.atmosphere.pressure_levels[end] * 
        buf.scene.atmosphere.pressure_unit
    
    d["SZA"] = buf.scene.solar_zenith * u"°"
    
    chi2 = G3RT.calculate_chi2(s)
    xgases = G3RT.calculate_xgas(buf.scene.atmosphere)

    # We only grab the first one, since we only have a single-band set-up
    d["chi2"] = first(values(chi2))
    
    # Get the XGAS values - remember these values have units
    for (k,v) in xgases
        d["X$(k)"] = v
    end

    # Get the column-averaged gravity
    d["g_mean"] = mean(
        buf.scene.atmosphere.gravity_levels 
        * buf.scene.atmosphere.gravity_unit
        )

    # Get the number of molecules for the total column (molec/m2)
    for atm in buf.scene.atmosphere.atm_elements
        if atm isa G3RT.GasAbsorber

            # = sum(num_dry_air_molecules * gas_vmr)
            d["$(atm.gas_name)"] = sum(
                first(values(buf.rt)).optical_properties.nair_dry
                .* G3RT.levels_to_layers((atm.vmr_levels * atm.vmr_unit) .|> NoUnits)
            ) * u"m^-2"

        end
    end
    
    d["ndry"] = sum(first(values(buf.rt)).optical_properties.nair_dry) * u"m^-2"

    d["iterations"] = G3RT.get_iteration_count(s)
    d["datetime"] = buf.scene.time
    d["unixtime"] = datetime2unix(d["datetime"])

    # Calculate RRMS (relative root-mean-square of residual)
    meas = G3RT.get_measured(s, swin);
    conv = buf.rt_buf.radiance.I[buf.rt_buf.indices[swin]];
    rrms = sqrt(mean((100 .* (conv .- meas) ./ maximum(meas)).^2))
    d["rrms"] = rrms * Unitful.percent

    # Pressure weighting "function"
    d["pressure_weights"] = G3RT.create_pressure_weights(buf.scene.atmosphere);

end



"""
Converts the big result dictionary into a DataFrame for easy storage
"""
function convert_to_df(
    d::Dict
)

    #=
        The big dictionary has follwoing structure:
            d[2]["./test_file.BIN"]["iteration_count"] = 3

        The first key corresponds to the window configuration index, 
        the second one corresponds to the file name of the processed spectrum,
        and the final key represents some result value. The second and third keys are
        dependent on the first one, meaning that there could be files processed for some
        first key, that weren't processed for the second key, and some results might be
        available for some configuration that are not available for another.

    =#


    
    # Since we are also going through most of the dictonary here,
    # might as well use the opportunity to store all result keys so
    # we know what data the final dataframe will hold.

    all_dfs = Dict{Integer, DataFrame}()

    for idx in keys(d)
        # Loop through the retrieval configurations


        # Find the data keys and their types
        this_d = first(values(d[idx]))
        # Extract the types of this directory
        type_dict = Dict(k => typeof(v) for (k,v) in this_d)
        unit_dict = Dict(k => v isa Unitful.Quantity ? unit(v) : 1 for (k,v) in this_d)

        # TODO
        # We should have some sensible ordering function here, let the datetime be first
        # xgases last, something along those lines..


        # Create an empty data frame with the right types
        this_df = DataFrame(
            [Vector{x[2]}() for x in type_dict],
            [Symbol(x[1]) for x in type_dict]
            )

        # Slot it into the bigger dictonary
        all_dfs[idx] = this_df

        #=
            Walk through each spectrum and push the results into
            the data frame
        =#

        for (fname, res) in d[idx]
            this_row = [d[idx][fname][n] for n in names(this_df)]
            push!(this_df, this_row)
        end

        # Sort the data frame by date
        sort!(this_df, "datetime")


        # Make unit replacements
        for (k, this_unit) in unit_dict

            if this_unit isa Unitful.Units
                transform!(this_df, k => ustrip => k)
                rename!(this_df, k => "$(k) [$(this_unit)]")
            end

        end
    end

    return all_dfs


end