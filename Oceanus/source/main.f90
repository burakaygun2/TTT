program main
    use mod_parameters
    use mod_variables
    use mod_preprocess
    use mod_preprocess_lsingle
    use mod_time_integration
    !! delete later
    use mod_SHTns
    use mod_matrix
    use mod_rhs
    use mod_solver
    use mod_coriolis
    implicit none
    ! integer :: j, m, jr, mr, irr
    ! integer :: il, ir, ilr, jm_shtns, info
    integer :: j
    print*, 'Starting Oceanus...'
    call initial_read
    print*, 'Calling preprocess...'
    call prep
    call write_tidal_potential
    print*, '#Preprocess is done...'
    if(SWITCH_vgradv == 1) call SHTns_init_o
    
    print*, 'Generating matrix with cut-off degree', max_jmax
    
    do j = 1, max_jmax
        call create_matrix(j)
    end do
    
    print*, forcing_k(jm(2,0))
    print*, forcing_k(jm(2,1))
    print*, forcing_k(jm(2,2))
    print*, 'compute scalar harm = ', compute_scalar_harm
    ! stop
    
    call time_integration_degree_2
    if(SWITCH_coriolis == 0) then
        call time_integration_static
    elseif((SWITCH_coriolis == 1).and.(SWITCH_vgradv == 0)) then
        call time_integration_coriolis
    elseif((SWITCH_coriolis == 1).and.(SWITCH_vgradv == 1)) then
        call time_integration_vgradv
    end if
    
end program main