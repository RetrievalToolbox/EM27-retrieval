function process_spectrum(
    meas, # Measured radiance dictionary (dispersion -> vector)
    time, # Time at measurement
    SZA, # Solar zenith angle
    psurf, # Surface pressure (with units)
    dispersion_coeffs, # Initial dispersion coefficients
    dispersions, # Dispersion dictionary (spectral window -> disperison)
    noise, # Noise sigma
    map_atm, # MAP atmosphere dictionary
    state_vector, # State vector
    isrf_dict, # ISRF dictionary
    buf; # Buffer
    max_iter=10
    )
   
    N_sv = length(state_vector)

    # Copy the contents of the MAP dictionary into the atmosphere, such
    # that we start we a fresh atmosphere every time.
    copy_map_to_atmosphere!(buf.scene, map_atm);

    # Create the retrieval grid, based on the surface pressure
    plevels = RE.create_UoL_pressure_grid(psurf, 400.0u"hPa", N_total=buf.scene.atmosphere.N_level)
    RE.ingest!(buf.scene.atmosphere, :pressure_levels, plevels)
        
    # Calculate altitude and gravity levels - this needs to be done internally, since
    # we must have the values for the finer MET grid, rather than the retrieval grid
    RE.calculate_altitude_and_gravity!(buf.scene)

    # Calculate layer quantities (from level quantities)
    RE.calculate_layers!(buf.scene.atmosphere)
    
    rt_buf = buf.rt_buf
    inst_buf = buf.inst_buf
    
    buf.scene.time = time
    buf.scene.solar_zenith = SZA
    buf.scene.solar_azimuth = 0.0

    # Loop through gases and set to MAP atmosphere values
    for atm in buf.scene.atmosphere.atm_elements
        
        if atm isa RE.GasAbsorber

            if atm.gas_name == "O2"
                # Make the O2 gas VMR vector
                o2_vmr = RE.atmospheric_profile_interpolator_linear(
                    ustrip.(Ref(buf.scene.atmosphere.pressure_unit), map_atm["Pressure"]),
                    ustrip.(Ref(atm.vmr_unit), map_atm["o2"]),
                    buf.scene.atmosphere.pressure_levels
                    )
            
                @views atm.vmr_levels[:] = o2_vmr 
            end

            if atm.gas_name == "H2O"
                # Make the H2O gas VMR vector
                h2o_vmr = RE.atmospheric_profile_interpolator_linear(
                    ustrip.(Ref(buf.scene.atmosphere.pressure_unit), map_atm["Pressure"]),
                    ustrip.(Ref(atm.vmr_unit), map_atm["h2o"]),
                    buf.scene.atmosphere.pressure_levels
                    )
            
                @views atm.vmr_levels[:] = h2o_vmr 
            end

            if atm.gas_name == "CO2"
                # Make the CO2 gas VMR vector (in ppm)
                co2_vmr = RE.atmospheric_profile_interpolator_linear(
                    ustrip.(Ref(buf.scene.atmosphere.pressure_unit), map_atm["Pressure"]),
                    ustrip.(Ref(atm.vmr_unit), map_atm["co2"]),
                    buf.scene.atmosphere.pressure_levels
                    )
            
                @views atm.vmr_levels[:] = co2_vmr 
            end

            if atm.gas_name == "CH4"
                # Make the CH4 gas VMR vector (in ppb)
                ch4_vmr = RE.atmospheric_profile_interpolator_linear(
                    ustrip.(Ref(buf.scene.atmosphere.pressure_unit), map_atm["Pressure"]),
                    ustrip.(Ref(atm.vmr_unit), map_atm["ch4"]),
                    buf.scene.atmosphere.pressure_levels
                    )
            
                @views atm.vmr_levels[:] = ch4_vmr 
            end            

            if atm.gas_name == "CO"
                # Make the CO gas VMR vector (in ppb)
                co_vmr = RE.atmospheric_profile_interpolator_linear(
                    ustrip.(Ref(buf.scene.atmosphere.pressure_unit), map_atm["Pressure"]),
                    ustrip.(Ref(atm.vmr_unit), map_atm["co"]),
                    buf.scene.atmosphere.pressure_levels
                    )
            
                @views atm.vmr_levels[:] = co_vmr 
            end 
            
        end
    end

    
    
    ##################################################
    # Prepare the SV with priors and prior covariances
    ##################################################

    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.SurfacePressureSVE)

        sve.first_guess = (
            buf.scene.atmosphere.pressure_levels[end] * 
            buf.scene.atmosphere.pressure_unit |> sve.unit) |> ustrip
        sve.prior_value = sve.first_guess
        sve.iterations[1] = sve.first_guess
        
    end
    

    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.SolarScalerPolynomialSVE)
        
        swin = sve.swin
        disp = dispersions[swin]
        rt = buf.rt[swin]
        this_meas = meas[disp]
        
        if sve.coefficient_order == 0

            # Need to reference the measured spectrum via the dispersion
            solar_idx = searchsortedfirst.(Ref(rt.solar_model.ww), (swin.ww_min, swin.ww_max))
            
            sve.first_guess = mean(this_meas[isfinite.(this_meas)]) / mean(rt.solar_model.irradiance[solar_idx[1]:solar_idx[2]])
            sve.prior_value = sve.first_guess
            sve.prior_covariance = (0.1 * sve.first_guess)^2

        end
    end
    
    # Ingest dispersion polynomial coefficient priors
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.DispersionPolynomialSVE)

        o = sve.coefficient_order
        if o < length(dispersion_coeffs)
    
            sve.first_guess = dispersion_coeffs[o+1]
            sve.prior_value = sve.first_guess

            if sve.coefficient_order == 0
                sve.prior_covariance = 1.0e-2
            else
                sve.prior_covariance = max(1e-5, (0.01 * sve.prior_value)^2)
            end
            
        end

    end
    
    # Empty out all state vector element iterations,
    # and populate with first guesses
    for sve in state_vector.state_vector_elements
        empty!(sve.iterations)
        push!(sve.iterations, sve.first_guess)
    end

    # Define the forward model here
    # (note that this is a simple "rebinding" of arguments,
    #  and is practically free performance-wise)
    fm!(SV) = forward_model!(
        SV,
        buf,
        inst_buf,
        rt_buf,
        dispersions,
        isrf_dict,
        0.0, #solar_doppler_factor,
        1.0, #solar_distance
    )

    # Create a view to the Sa buffer, such that the
    # solver can use it later on.
    Sa = Diagonal(RE.get_prior_covariance(state_vector))

    # Calculate measurement noise dictionary
    noise_dict = Dict{RE.AbstractDispersion, AbstractVector}()
    for disp in values(dispersions)
        meas_noise = similar(meas[disp])
        @views meas_noise[:] .= noise
        noise_dict[disp] = meas_noise
    end
       
    # Create a solver object for our chosen solver
    solver = RE.IMAPSolver(
        fm!,
        state_vector,
        Sa,
        max_iter,
        1.0,
        dispersions,
        rt_buf.indices,
        rt_buf.radiance,
        rt_buf.jacobians,
        meas,
        noise_dict,
    )

    # Perform the iterations
    iter_count = 0
    converged = false
    while (iter_count < solver.max_iterations)

        converged = RE.check_convergence(solver)
        
        if converged
            @debug "Reached convergence!"
            break
        end

        iter_success = false

        try
            iter_success = RE.next_iteration!(solver)
            iter_count += 1
        catch e
            @error "FAILED executing iteration with this error:"
            @error e
            break
        end

        if !iter_success
            @error "FAILED iteration #$(iter_count)"
            break
        else
            @debug "Successful iteration #$(iter_count)"
        end

        chi2 = RE.calculate_chi2(solver)

        msg = "#$(iter_count): $(chi2)"
        #@info msg

    end
    @debug "Iterations done."

    # Update atmosphere one last time after the final state vector update was done
    RE.atmosphere_element_statevector_update!(buf.scene.atmosphere.atm_elements, state_vector)
    RE.atmosphere_statevector_update!(buf.scene.atmosphere, state_vector)
    RE.update_specific_humidity_from_H2O!(buf.scene.atmosphere)

    #RE.calculate_OE_quantities(solver)
    #chi2 = RE.calculate_chi2(solver)

    #msg = @sprintf("%d, it#: %d, Chi2: %.4f, converged: %s", snid, iter_count, chi2, converged)
    #@debug msg

    return solver
    
end
