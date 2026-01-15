function forward_model!(
    SV::RE.AbstractStateVector,
    buf::RE.EarthAtmosphereBuffer,
    inst_buf::RE.InstrumentBuffer,
    rt_buf::RE.ScalarRTBuffer,
    dispersions::Dict{<:RE.AbstractSpectralWindow, <:RE.AbstractDispersion},
    isrf_dict::Dict{<:RE.AbstractDispersion, <:RE.AbstractISRF},
    solar_doppler_factor::AbstractFloat,
    solar_distance::AbstractFloat, # in multiples of AU
    )


    # Update the dispersions object based on the state vector
    # (You must make sure that the dispersion state vector elements
    #  are bound to the used dispersion!)
    for disp in values(dispersions)
        RE.update_dispersion!(disp, SV)
    end
    
    # Re-calculate indices for all windows
    RE.calculate_indices!(buf)

    # Ingest retrieved H2O into specific humidity!
    RE.update_specific_humidity_from_H2O!(buf.scene.atmosphere)

    # Update the atmosphere with SV contents (e.g. temperature scale)
    RE.atmosphere_statevector_update!(buf.scene.atmosphere, SV)

    for (idx, sve) in RE.StateVectorIterator(SV, RE.SurfacePressureSVE)

        new_psurf = ustrip(
            buf.scene.atmosphere.pressure_unit,
            RE.get_current_value_with_unit(sve)
            )
        
        @info "Updating PSURF to $(new_psurf)"
        buf.scene.atmosphere.pressure_levels[end] = new_psurf

        # .. and calculate the mid-layer values that need to be
        # inserted into the atmosphere object.
        @views buf.scene.atmosphere.pressure_layers[:] =
            RE.levels_to_layers(buf.scene.atmosphere.pressure_levels)

        # Re-calculate z and g
        RE.calculate_altitude_and_gravity!(buf.scene)

    end
    
    # Update the solar scaler field in the RT container
    for rt in values(buf.rt)
        RE.solar_scaler_statevector_update!(rt, SV)
    end

    # Obtain down-sampled solar radiance at the required wavelength grid
    for swin in keys(buf.rt)

        rt = buf.rt[swin]
        
        RE.calculate_solar_irradiance!(
            rt,
            swin,
            rt.solar_model,
            doppler_factor=solar_doppler_factor
        )

        # Scale solar irradiance according to relative solar distance
        @views rt.hires_solar.I[:] /= (solar_distance^2)
        
    end


    
    # Given the current statevector - modify any atmospheric element
    # according to whatever the state vector element dictates.
    # (e.g. scale the VMR of a gas according to the SV)
    RE.atmosphere_element_statevector_update!(buf.scene.atmosphere.atm_elements, SV)
    
    # Calculate optical properties inside buffer
    RE.calculate_earth_optical_properties!(buf, SV, N_sublayer=1)
    
    # Perform RT
    for rt in values(buf.rt)
        RE.calculate_radiances_and_jacobians!(
            rt,
            buf.scene.observer
        )
    end
    
    # Apply the ISRF to the hi-res spectrum to obtain
    # detector-level radiance
    for swin in keys(buf.rt)
        
        disp = dispersions[swin]

        success = RE.apply_isrf_to_spectrum!(
                inst_buf, # Convolution buffer
                isrf_dict[disp], # ISRF table
                disp, # Dispersion
                buf.rt[swin].hires_radiance.I, # radiance
                swin
            )

        if !success
            @debug "Application of ISRF on radiances FAILED for window $(swin)."
            return false
        else
            @debug "Application of ISRF on radiances successful for window $(swin)."
        end

        # Store in RT buffer, but account for unit differences!
        @views rt_buf.radiance.I[rt_buf.indices[swin]] = inst_buf.low_res_output[disp.index] * 
            buf.rt[swin].radiance_unit / rt_buf.radiance_unit 
            
    end


    # Apply any radiance correction, such as ZLO
    for sve in SV.state_vector_elements
        RE.apply_radiance_correction!(rt_buf, sve)
    end

    # Apply the ISRF to the hires Jacobians
    for (i, sve) in enumerate(SV.state_vector_elements)
        
        if RE.calculate_jacobian_before_isrf(sve)

            for swin in keys(dispersions)

                disp = dispersions[swin]
                
                success = RE.apply_isrf_to_spectrum!(
                    inst_buf,
                    isrf_dict[disp],
                    disp,
                    buf.rt[swin].hires_jacobians[sve].I,
                    swin
                )
    
                if !success
                    @debug "Application of ISRF on Jacobian $(i) for $(sve) FAILED."
                    return false
                else
                    @debug "Application of ISRF on Jacobian $(i) for $(sve) successful."
                end
    
                @views rt_buf.jacobians[sve].I[rt_buf.indices[swin]] = inst_buf.low_res_output[disp.index]
            end
        end
    end

    # Post-ISRF Jacobian calculations

    # ZLO jacobians
    for (i, sve) in RE.StateVectorIterator(SV, RE.ZeroLevelOffsetPolynomialSVE)
        RE.calculate_jacobian!(rt_buf, sve)          
    end
    
    for (i, sve) in RE.StateVectorIterator(SV, RE.DispersionPolynomialSVE)

        for swin in keys(dispersions)

            disp = dispersions[swin]
            
            success = RE.calculate_dispersion_polynomial_jacobian!(
                    inst_buf,
                    sve,
                    isrf_dict[disp],
                    buf.rt[swin].hires_radiance.I,
                )

            rt_buf.jacobians[sve].I[rt_buf.indices[swin]] = inst_buf.low_res_output[disp.index]
        
            if !success
                @debug "Calculation of Jacobian for $(sve) FAILED."
                return false
            else
                @debug "Calculation of Jacobian for $(sve) successful."
            end

        end

    end
    

    #for sve in keys(rt_buf.jacobians)
    #    @views rt_buf.jacobians[sve].I[rt_buf.indices[swin]] *= 1.0 * buf.rt[swin].radiance_unit / rt_buf.radiance_unit
    #end

    # Important!
    # Re-set the atmosphere to its prior state!
    RE.atmosphere_element_statevector_rollback!(buf.scene.atmosphere.atm_elements, SV)
    RE.atmosphere_statevector_rollback!(buf.scene.atmosphere, SV)

    return true

end