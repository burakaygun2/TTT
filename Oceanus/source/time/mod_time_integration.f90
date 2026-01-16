module mod_time_integration
    use mod_parameters
    use mod_variables
    use mod_allocation
    use mod_displacement
    use mod_matrix
    use mod_solver
    use mod_potential
    use mod_norms
    use mod_coriolis
    use mod_write_f
    use mod_write_h5
    use mod_SHTns
    use mod_math
    use mod_preprocess
    implicit none
    logical :: WRITE_velocity_cond
    real(dp) :: ocean_vol, reynolds, rossby, ocean_thick, ocean_dens
contains

subroutine time_integration_degree_2
    implicit none
    integer :: it
    real(dp) :: bt, ft
    ocean_vol = Pi * (4.0_dp/3.0_dp) * (layer(1)%r_i(layer(1)%n)**3 - layer(1)%r_i(1)**3)
    ocean_thick = layer(1)%r_i(layer(1)%n) - layer(1)%r_i(1)
    ocean_dens = 0.5_dp * (layer(1)%rho_i(1)+layer(1)%rho_i(layer(1)%n))
    call overwritejmax(2)
    t = time_step; it = 1
    love_k = -1.0_dp
    bt  = omp_get_wtime()
    do  
        call solver( t)
        call boundary_displacement
        call ind_pot2(it)
        call self_gravity
        call velocity_norm
        reynolds = ocean_dens * ocean_thick * sqrt(layer(1)%v_norm/ocean_vol) / layer(1)%eta_i(1)
        rossby   = sqrt(layer(1)%v_norm/ocean_vol) / (2.0_dp * ang_vel * ocean_thick)
        ! print*, reynolds, rossby, layer(1)%eta_i(1) / (2.0_dp*ocean_dens*ang_vel*(ocean_thick)**2), sqrt(layer(1)%v_norm/ocean_vol)
        ! call determine_dt
        ! call find_boundary_layer_thickness(it)
        it = it + 1; t = t + time_step
        if (it == 10*num_time_steps+1) exit
    end do
    ft  = omp_get_wtime()
    print*, '=============================================================' 
    print'(A, F20.8, A)', 'Initial runtime: ', ft - bt, 's'
    print*, '=============================================================' 
end subroutine time_integration_degree_2


subroutine time_integration_coriolis
    implicit none
    integer :: it
    real(dp) :: bt, ft

    call setup_time_array
    call open_file_full_potential
    call overwritejmax(init_jmax)
    t = time_step; it = 1
    bt  = omp_get_wtime()
    do  
        call solver( t)
        call coriolis3_r
        call boundary_displacement
        call ind_pot2_coriolis(it)
        call self_gravity
        call velocity_norm
        call compute_dissipation(it)
        if (it>= (final_period-2)*num_time_steps) call write_file_full_potential(t)
        if(mod(it, num_time_steps) == 0) call period_update(it)
        call eval_write_cond(it)
        it = it + 1; t = t + time_step
        if (it == nint(final_period*num_time_steps)+1) exit
    end do
    ft  = omp_get_wtime()
    print*, '=============================================================' 
    print'(A, F15.8, A)', 'TOTAL time: ', ft - bt, 's'
    print*, '=============================================================' 
end subroutine time_integration_coriolis

subroutine time_integration_static
    implicit none
    integer :: it
    real(dp) :: bt, ft

    call setup_time_array
    call overwritejmax(2)
    t = time_step; it = 1
    bt  = omp_get_wtime()
    
    do  
        call solver( t)
        call boundary_displacement
        call ind_pot2(it)
        call self_gravity
        call velocity_norm
        call compute_dissipation(it)
        call write_file_full_potential(t)
        if(mod(it, num_time_steps) == 0) call period_update(it)
        call eval_write_cond(it)
        it = it + 1; t = t + time_step
        if (it == nint(final_period*num_time_steps)+1) exit
    end do
    ft  = omp_get_wtime()
    print*, '=============================================================' 
    print'(A, F15.8, A)', 'TOTAL time: ', ft - bt, 's'
    print*, '=============================================================' 
end subroutine time_integration_static

subroutine time_integration_vgradv
    implicit none
    integer :: it
    real(dp) :: bt, ft

    call setup_time_array
    call overwritejmax(init_jmax)
    t = time_step; it = 1
    bt  = omp_get_wtime()
    open(600, file='test.dissipation')
    do  
        call solver(t)
        call coriolis3_r
        call vgradv
        call boundary_displacement
        call ind_pot2_coriolis(it)
        call self_gravity
        call velocity_norm
        call compute_dissipation(it)
        if(it > 28*num_time_steps) then
            if(mod(it, 5) == 0) print*, layer(1)%spectrum_v(:)
        end if
        write(600, *) t, layer(1)%dissipation*(2*layer(1)%eta_i(1)), layer(1)%v_norm
        if(mod(it, num_time_steps) == 0) then
            print*, 'Reached one period at ', t 
            call period_update(it)
        end if
        call eval_write_cond(it)
        it = it + 1; t = t + time_step
        if (it == nint(final_period*num_time_steps)+1) exit
    end do
    ft  = omp_get_wtime()
    print*, '=============================================================' 
    print'(A, F15.8, A)', 'Time stepping: ', ft - bt, 's'
    print*, '=============================================================' 
    close(600)
    
end subroutine time_integration_vgradv

subroutine period_update(it2)
    implicit none
    integer :: ip, it2
    ip = it2/num_time_steps
    call write_time_series_h5(it2)
    call write_periodic_h5(ip)
    time = time + rotation_period
end subroutine period_update

subroutine eval_write_cond(it2)
    implicit none
    integer, intent(in) :: it2
    WRITE_velocity_cond = ((WRITE_velocity == 1).and.(it2/num_time_steps >= start_output).and.(mod(it2, write_output)==0))
    if (WRITE_velocity_cond) call write_velocity_shc_h5(it2)
end subroutine


end module mod_time_integration