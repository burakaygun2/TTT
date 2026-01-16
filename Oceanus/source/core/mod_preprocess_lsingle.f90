module mod_preprocess_lsingle
    use mod_parameters
    use mod_variables
    use mod_allocation
    use mod_read_f
    use mod_write_f
    use mod_write_h5
    use mod_preprocess
    implicit none
    
contains
    
    subroutine prep_lsingle
        implicit none
        
        call prepare_single_layer
        call single_complete_nr
        call assign_lsingle
        if(SWITCH_nondim == 'n') call single_setup_nondim
        call mid_layers_single
        call make_dir
        call write_summary_lsingle_h5
        call write_radial_profiles_lsingle_h5
        call generate_timeseries_space_h5
        call generate_velocity_shc_h5
        call alloc_NM
        call alloc_RHS
        call alloc_vars_lsingle
        print*, 'Single layer is prepared'
    end subroutine prep_lsingle

    subroutine prepare_single_layer
        implicit none
        lsingle%n = N_initial
        call alloc_single_layer
        call alloc_rads_lsingle
        call complete_time
        print*, 'Non-dimensional time step...', time_step

    end subroutine prepare_single_layer

    subroutine single_setup_nondim
        implicit none
        lsingle%g_i   = 1.0_dp
        lsingle%rho_i = 1.0_dp
        lsingle%eta_i = 10.0_dp ** (Ocean_viscosity) !! Ekman number
        ang_vel       = 1.0_dp
        lsingle%r_i   = lsingle%r_i / maxval(lsingle%r_i)
    end subroutine single_setup_nondim


    subroutine single_complete_nr
        implicit none

        lsingle%sm = 7*lsingle%n + 2
        lsingle%tm = 3*lsingle%n + 1

        smtot      = lsingle%sm
        tmtot      = lsingle%tm

        total_scalar_harm = jm (Max_jmax, mmax)
        total_vector_harm = jml(Max_jmax, mmax, Max_jmax+1)

        compute_scalar_harm =  jm(init_jmax, mmax)
        compute_vector_harm = jml(init_jmax, mmax, init_jmax+1)

        call jmindx_
    end subroutine

    subroutine mid_layers_single
        implicit none
        
        call interface2layer(lsingle%r_i  , lsingle%r_l  )
        call interface2layer(lsingle%rho_i, lsingle%rho_l)
    end subroutine mid_layers_single

    subroutine assign_lsingle
        implicit none
        integer :: i, ir

        i = 1
        do ir = 1, lsingle%n
            lsingle%r_i(ir)   = temp_radius  (i)
            lsingle%rho_i(ir) = temp_density (i)
            i = i + 1
        end do
    end subroutine assign_lsingle
end module mod_preprocess_lsingle