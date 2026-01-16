module mod_preprocess
    use mod_parameters
    use mod_variables
    use mod_allocation
    use mod_read_f
    use mod_write_f
    use mod_write_h5
    implicit none
    
contains


subroutine prep
    implicit none
    
    ! call initial_read
    call number_of_layers
    call alloc_rads
    call assign_layers
    call find_ocean_layer
    call complete__
    call SWITCH_imp
    call eccentricity_time
    if(SWITCH_hydrostatic_core == 'h') call hydrostatic_core
    call mid_layers
    call radial_grads
    if(SWITCH_one_layer == 'm') then
        call compute_mass
        call compute_grav_acc
    else
        call compute_mass_grav_acc_single
    end if
    if(SWITCH_nondim == 'n') call set_nondim
    call compute_MOI
    call make_dir
    call write_summary_h5
    call write_radial_profiles_h5
    call generate_timeseries_space_h5
    call generate_velocity_shc_h5
    call alloc_NM
    call alloc_RHS
    call alloc_vars
    call write_initial_signal
    call read_initial_cond
    call determine_mindr
    call max_forcing
    call setup_time_array

end subroutine prep



subroutine initial_read
    implicit none

    call read_input
    call read_radial_profile
    if(unit_period == 'd') rotation_period = rotation_period * day2sec
    obl = obl * deg_to_rad
    if(temp_radius(1) - temp_radius(2)>0.0) call reverse_temps
    temp_shear_m = temp_shear_m * GPa_to_Pa
    temp_bulk_m  = temp_bulk_m  * GPa_to_Pa
end subroutine initial_read

subroutine number_of_layers
    implicit none
    integer :: il, ir, lnn(1:150), i
    real(dp) :: rr(1:150)
    logical :: cond1
    
    il=0; ir=1;
    do i = 1, N_initial-1
        cond1 = abs(temp_radius(i+1) - temp_radius(i)) < 1d-8
        if(cond1) then
            il = il + 1
            rr(il) = temp_radius(i)
            lnn(il) = ir
            ir = 0
        end if
        ir = ir + 1
    end do
    numoflayers = il + 1;            !if(numoflayers/=5) stop "Number of layers do not correspond to Ganymede"
    lnn(il+1)   = ir
    rr(il+1)    = temp_radius(N_initial)
    print*, numoflayers

    allocate(layer(numoflayers), rad(0:numoflayers), layer_ocean(numoflayers), layer_ocean_id(numoflayers))
    print*, 'Number of layers...', numoflayers
    rad(0) = 0.0_dp
    do il = 1, numoflayers
        rad(il)     = rr(il)
        layer(il)%n = lnn(il)
    end do

end subroutine number_of_layers

subroutine assign_layers
    implicit none
    integer :: i, il, ir, nn

    i = 1
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn
            layer(il)%r_i(ir)   = temp_radius  (i)
            layer(il)%rho_i(ir) = temp_density (i)
            layer(il)%K_i(ir)   = temp_bulk_m  (i)
            layer(il)%mu_i(ir)  = temp_shear_m (i)
            layer(il)%eta_i(ir) = temp_visc    (i)
            i = i + 1
        end do
    end do
end subroutine assign_layers



subroutine complete__
    implicit none
    call complete_time
    call complete_nr
end subroutine complete__

subroutine complete_time
    implicit none
    
    time_step = rotation_period / real(num_time_steps, dp)
    if (SWITCH_nondim == 'n') then
        ang_vel = rotation_period
    else
        ang_vel = 2.0_dp * Pi / rotation_period
    end if

    dtmin_r = 1e30_dp; dtmin_t = 1e30_dp
    change_time_step = 0
end subroutine complete_time

subroutine complete_nr
    implicit none
    integer :: il

    do il = 1, numoflayers
        layer(il)%sm = layer(il)%n * 7 + 2
        layer(il)%tm = layer(il)%n * 3 + 1
    end do

    Nlay = 0; Nrad = 0
    do il = 1, numoflayers
        Nlay = Nlay + layer(il)%n + 1
        Nrad = Nrad + layer(il)%n
    end do

    !!this is for velocity and Coriolis term
    layer(1)%lbound = 1
    layer(1)%ubound = layer(1)%n + 1

    do il = 2, numoflayers
        layer(il)%lbound = (il-1)*layer(il-1)%ubound + 1
        layer(il)%ubound = layer(il-1)%ubound + layer(il)%n + 1
    end do

    
    smtot = sum(layer(1:numoflayers)%sm)
    tmtot = sum(layer(1:numoflayers)%tm)

    total_scalar_harm = jm (Max_jmax,mmax)
    total_vector_harm = jml(Max_jmax,mmax,Max_jmax+1)
    compute_scalar_harm =  jm(init_jmax, mmax)
    compute_vector_harm = jml(init_jmax, mmax, init_jmax+1)

    ! if(SWITCH_vgradv == 1) then
    !     total_scalar_harm = jm (Max_jmax,Max_jmax)
    !     total_vector_harm = jml(Max_jmax,Max_jmax,Max_jmax+1)
    !     compute_scalar_harm =  jm(init_jmax, init_jmax)
    !     compute_vector_harm = jml(init_jmax, init_jmax, init_jmax+1)
    ! end if

    call jmindx_
    
    jstart = 2; jfinal = 2; status = 0
end subroutine complete_nr


subroutine reverse_temps
    implicit none
    call reverse_array(temp_radius  (:))
    call reverse_array(temp_density (:))
    call reverse_array(temp_bulk_m  (:))
    call reverse_array(temp_shear_m (:))
    call reverse_array(temp_visc    (:))
end subroutine reverse_temps


subroutine mid_layers
    implicit none
    integer :: il

    do il = 1, numoflayers
        call interface2layer(layer(il)%r_i  , layer(il)%r_l  )
        call interface2layer(layer(il)%rho_i, layer(il)%rho_l)
    end do
    if(SWITCH_one_layer == 'm')layer(1)%r_l(1) = 0.0_dp
end subroutine mid_layers

subroutine radial_grads
    implicit none
    integer :: il, ir, nn
    real(dp) :: dr, drho

    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn
            dr   = layer(il)%r_l  (ir+1) - layer(il)%r_l  (ir)
            drho = layer(il)%rho_l(ir+1) - layer(il)%rho_l(ir)
            layer(il)%rho_r_i(ir) = drho / dr
        end do
    end do

    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 2, nn
            dr   = layer(il)%r_i  (ir  ) - layer(il)%r_i  (ir-1)
            drho = layer(il)%rho_i(ir  ) - layer(il)%rho_i(ir-1)
            layer(il)%rho_r_i(ir) = drho / dr
        end do
    end do
end subroutine radial_grads

subroutine check_layers
    implicit none
    integer  :: il

    do il = 1, numoflayers
        if (maxval(abs(layer(il)%rho_r_i))<1e-8_dp) then
            call delloc_real1(layer(il)%rho_r_i)
            call delloc_real1(layer(il)%rho_r_l)
            constant_density = .TRUE.
        end if
    end do
end subroutine check_layers

subroutine SWITCH_imp
    implicit none
    integer :: il, oceanindex


    !! The rocky parts are always incompressible
    if ((SWITCH_compressible == 1).and.(SWITCH_one_layer == 'm')) then
        layer(1)%K_i = 1e30_dp
        layer(2)%K_i = 1e30_dp
    end if

    if(SWITCH_ocean_status == 'h') then
        do il = 1, numoflayers
            if(layer_ocean_id(il) == 1) then
                layer(il)%mu_i  = 1e-2_dp
                layer(il)%eta_i = 1e30_dp
            end if
        end do
    else
        do il = 1, numoflayers
            do oceanindex = 1, num_oceans
                if(il == layer_ocean(oceanindex)) then
                    print*, 'ocean index is ', il
                    layer(il)%mu_i  = 1e30_dp
                    layer(il)%eta_i = 10.0_dp**layer(il)%eta_i
                    if(SWITCH_override_viscosity == 1) then
                        layer(il)%eta_i = 10.0_dp**ocean_viscosity
                    end if
                end if
            end do
        end do
    end if

    if(SWITCH_solid_layer == 'e') then
        do il = 1, numoflayers
            do oceanindex = 1, num_oceans
                if(il /= layer_ocean(oceanindex)) layer(il)%eta_i = 1e30_dp
            end do
        end do
    else
        do il = 1, numoflayers
            if(layer_ocean_id(il) /= 1) then
                layer(il)%eta_i = 10.0_dp ** layer(il)%eta_i
            end if
        end do
    end if
end subroutine SWITCH_imp

subroutine find_ocean_layer
    implicit none
    integer :: il, oceanindex, ils
    layer_ocean_id = 0
    layer_ocean = 0
    oceanindex = 1
    if ((SWITCH_hydrostatic_core == 'h').and.(SWITCH_one_layer == 'm')) then
        ils = 2
    else
        ils = 1
    end if
    do il = ils, numoflayers
        if(((minval(layer(il)%eta_i)) < 20.0_dp).and.((maxval(layer(il)%mu_i)) < 10.0_dp)) then
            layer_ocean_id(il) = 1
            layer_ocean(oceanindex) = il
            oceanindex = oceanindex + 1
        end if
    end do
    num_oceans = oceanindex - 1
    print*, 'Number of viscous liquid layers...', oceanindex-1
end subroutine find_ocean_layer

subroutine hydrostatic_core
    implicit none
    print*, 'hydrostatic core is set.'
    layer(1)%mu_i  = 1e2_dp
    layer(1)%eta_i = 1e30_dp
end subroutine hydrostatic_core

subroutine compute_mass
    implicit none
    integer :: il, ir, nn
    real(dp) :: f, vol

    f   = 4.0_dp * Pi / 3.0_dp
    vol = f * layer(1)%r_i(1)
    layer(1)%mass_i(1) = layer(1)%rho_i(1)*vol
    do il = 1, numoflayers
        nn = layer(il)%n
        if(il >= 2) layer(il)%mass_i(1) = layer(il-1)%mass_i(layer(il-1)%n)
        do ir = 2, layer(il)%n
            vol = f * (layer(il)%r_i(ir)**3 - layer(il)%r_i(ir-1)**3)
            layer(il)%mass_i(ir) = layer(il)%mass_i(ir-1) + vol *layer(il)%rho_i(ir)
        end do
        call interface2layer(layer(il)%mass_i, layer(il)%mass_l)
    end do
    layer(1)%mass_l(1) = 0.0_dp
end subroutine compute_mass

subroutine compute_grav_acc
    implicit none
    integer :: il,ir, sp

    do il = 1, numoflayers
        do ir = 1, layer(il)%n
            layer(il)%g_i(ir) = layer(il)%mass_i(ir) * Grav_Const / layer(il)%r_i(ir)**2
        end do
    end do
    do il = 1, numoflayers
        if(il==1)sp=2
        if(il/=1)sp=1
        do ir = sp, layer(il)%n+1
            layer(il)%g_l(ir) = layer(il)%mass_l(ir) * Grav_Const / layer(il)%r_l(ir)**2
        end do
    end do
    layer(1)%g_l(1) = 0.0_dp
end subroutine compute_grav_acc

subroutine compute_mass_grav_acc_single
    implicit none
    integer  :: ir
    real(dp) :: f, vol

    f   = 4.0_dp * Pi / 3.0_dp
    vol = f * layer(1)%r_i(1)
    layer(1)%mass_i(1) = bottom_grav_acc * layer(1)%r_i(1)**2 / Grav_Const
    layer(1)%g_i(1)    = bottom_grav_acc
    do ir = 2, layer(1)%n
        vol = f * (layer(1)%r_i(ir)**3 - layer(1)%r_i(ir-1)**3)
        layer(1)%mass_i(ir) = layer(1)%mass_i(ir-1) + vol *layer(1)%rho_i(ir)
        layer(1)%g_i(ir) = layer(1)%mass_i(ir) * Grav_Const / layer(1)%r_i(ir)**2
    end do
end subroutine compute_mass_grav_acc_single

subroutine compute_MOI
    implicit none
    integer :: i, ir
    real(dp) :: fac, dr
    
    fac=8.0_dp*pi/(3.0_dp*(rad(numoflayers)**2)*layer(numoflayers)%mass_i(layer(numoflayers)%n))

    moi = 0.0_dp
    do i = 1, numoflayers
        do ir = 1, layer(i)%n-1
            dr = layer(i)%r_i(ir+1) - layer(i)%r_i(ir)
            moi = moi + (layer(i)%rho_i(ir+1)*layer(i)%r_i(ir+1)**4 + layer(i)%rho_i(ir)*layer(i)%r_i(ir)**4)*0.5_dp*dr
        end do
    end do
    moi = moi*fac
end subroutine compute_MOI

subroutine determine_mindr
    implicit none
    integer :: ir, il, nn

    min_dr = 1d20
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn
            if ( abs(layer(il)%r_l(ir+1) - layer(il)%r_l(ir)) < min_dr) then
                min_dr = abs(layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
            end if
        end do
    end do
end subroutine determine_mindr

subroutine max_forcing
    implicit none
    integer  :: it
    real(dp) :: tt, rr, rk, ik, rk1, ik1
    complex(dp) :: a, b
    rk = 0.0; ik = 0.0; rk1 = 0.0; ik1 = 0.0
    if(SWITCH_mmtides == 2) then
        forcing_k = dzero
        rr = layer(numoflayers)%r_i(layer(numoflayers)%n)
        do it = 1, num_time_steps
            tt = (it-1)*time_step
            a = lunar_tides(1.0_dp, rr, 0, tt, 'p')
            if(real(a) > real(forcing_k(jm(2,0)))) forcing_k(jm(2,0)) = a
            a = lunar_tides(1.0_dp, rr, 2, tt, 'p')
            if(real(a) > rk) rk = real(a)
            if(imag(a) > ik) ik = imag(a)
        end do
        forcing_k(jm(2,2)) = cmplx(rk, ik, dp)
    elseif(SWITCH_mmtides == 0) then
        forcing_k = dzero
        rr = layer(numoflayers)%r_i(layer(numoflayers)%n)
        do it = 1, num_time_steps
            tt = (it-1)*time_step
            a = tides(1.0_dp, rr, 0, tt, 'p')
            if(real(a) > real(forcing_k(jm(2,0)))) forcing_k(jm(2,0)) = a
            a = tides(1.0_dp, rr, 2, tt, 'p')
            if(real(a) > rk) rk = real(a)
            if(imag(a) > ik) ik = imag(a)
            b = tides(1.0_dp, rr, 1, tt, 'p')
            if(real(b) > rk1) rk1 = real(b)
            if(imag(b) > ik1) ik1 = imag(b)
        end do
        forcing_k(jm(2,2)) = cmplx(rk , ik , dp)
        forcing_k(jm(2,1)) = cmplx(rk1, rk1, dp)
    end if
end subroutine max_forcing


subroutine jmindx_
    implicit none
    integer :: j, m, l

    allocate(jmindx(total_scalar_harm, 3), jmlindx(total_vector_harm, 4))
    do j = 0, max_jmax
        do m = 0, min(mmax, j)
            jmindx(jm(j,m), :)= [jm(j,m), j, m]
        end do
    end do

    do j = 0, Max_jmax
        do m = 0, min(mmax, j)
            do l = abs(j-1), j+1
                jmlindx(jml(j,m,l), :) = [jml(j,m,l), j, m, l]
            end do
        end do
    end do
end subroutine jmindx_

subroutine setup_time_array
    implicit none
    integer :: it

    do it = 1, num_time_steps
        time(it) = it*time_step
    end do
end subroutine setup_time_array


subroutine set_nondim
    implicit none
    layer(1)%g_i   = layer(1)%g_i/layer(1)%g_i(layer(1)%n)
    layer(1)%rho_i = 1.0_dp
    layer(1)%eta_i = 10.0_dp ** (Ocean_viscosity) !! Ekman number
    layer(1)%r_i   = layer(1)%r_i / maxval(layer(1)%r_i)
    call interface2layer(layer(1)%r_i  , layer(1)%r_l  )
    call interface2layer(layer(1)%rho_i, layer(1)%rho_l)
end subroutine set_nondim


subroutine eccentricity_time
    implicit none

    call read_eccentricity
end subroutine eccentricity_time
end module mod_preprocess