module mod_write_f
    use mod_parameters
    use mod_variables
    use mod_potential
    use hdf5
    implicit none
contains

subroutine make_dir
    implicit none
    
    call system("mkdir -p "//trim(dir_output))
    call system("mkdir -p "//trim(dir_output)//"/sh_coeffs/density")
    call system("mkdir -p "//trim(dir_output)//"/sh_coeffs/velocity")
    call system("mkdir -p "//trim(dir_output)//"/sh_coeffs/stress")
    call system("mkdir -p "//trim(dir_output)//"/sh_coeffs/displacement")
    call system("mkdir -p "//trim(dir_output)//"/on_grid/density")
    call system("mkdir -p "//trim(dir_output)//"/on_grid/velocity")
    call system("mkdir -p "//trim(dir_output)//"/on_grid/stress")
    call system("mkdir -p "//trim(dir_output)//"/on_grid/displacement")
    call system("mkdir -p "//trim(dir_output)//"/checkpoint")
end subroutine make_dir


subroutine write_summary
    implicit none
    integer :: i
    character(200) :: fmt, file_summary
    write(file_summary, '(A)') trim(dir_output)//'/summary.txt'
    open(200, file=file_summary)
    write(200,'(A)') "Material Parameters:aaaa "
    fmt = "(A20, A20, A20,A20, A20, A20, A20)"
    write(200,fmt) "Layer number", " Radius (km)", " Density (kg/m3)", " Viscosity (Pa s)", " Shear Mod. (Pa)", "Bulk Mod. (Pa)", " Grav. Acc. (m/s2)"
    fmt = "(I20, F20.3, F20.3, ES20.3E2, ES20.3E2, ES20.3E2, F20.3)"
    do i = 1, numoflayers
        write(200, fmt) &
                        i, rad(i)*1d-3        , &
                        layer(i)%rho_i(layer(i)%n)          , &
                        layer(i)%eta_i(layer(i)%n)          , &
                        layer(i)%mu_i(layer(i)%n)           , &
                        layer(i)%K_i (layer(i)%n)           , &
                        layer(i)%g_i(layer(i)%n)
    end do
    write(200, '(A)') "======================================================================"
    write(200, '(A36,ES10.3E2)') 'Total mass of the body (kg)       = ', layer(numoflayers)%mass_i(layer(numoflayers)%n)
    write(200, '(A36,ES10.3E2)') 'Gravitational acc.     (m/s^2)    = ', layer(numoflayers)%g_i(layer(numoflayers)%n)
    write(200, '(A36,F10.5)')    'MOI factor                        = ', moi
    write(200, '(A36,F10.3)')    'Thinkness of the Core (km)        = ', (rad(1))*m_to_km
    write(200, '(A36,F10.3)')    'Thinkness of the  BMO (km)        = ', (rad(2)-rad(1))*m_to_km
    write(200, '(A36,F10.4)' )   'Eccentricity                      = ', ecc
    write(200, '(A36,F10.2)')    'Period (s)                        = ', rotation_period
    write(200, '(A36,ES10.3E2)') 'Angular velocity                  = ', ang_vel
    write(200, '(A36,ES10.3E2)') 'Core Ekman number                 = ', layer(1)%eta_i(1) / (2.0_dp*layer(1)%rho_i(1)*ang_vel*rad(1)**2)
    write(200, '(A36,ES10.3E2)') 'BMO  Ekman number                 = ', layer(2)%eta_i(1) / (2.0_dp*layer(2)%rho_i(1)*ang_vel*(rad(2)-rad(1))**2)
    write(200, '(A36,F10.2)')    'Mean density (kg/m3)              = ', layer(numoflayers)%mass_i(layer(numoflayers)%n)/(4.0_dp*pi*(rad(numoflayers)**3)/3.0_dp)
    write(200, '(A36,  A36)')    'Status of Ocean                   = ', trim(SWITCH_ocean_status)
    write(200, '(A36,  A36)')    'Status of Solid                   = ', trim(SWITCH_solid_layer)
    write(200, '(A)') '  '
    write(200, '(A)') "======================================================================"
    write(200, '(A)') "======================================================================"
    write(200, '(A)') '  '
    write(200, '(A)') "Numerical Parameters:"
    fmt = "(A20, A20, A20, A20)"
    write(200,fmt) "Layer number", "Number of layers", "Coriolis effect", "Inertia"
    fmt = "(I20, I20,  I20, I20)"
    do i = 1, numoflayers
        write(200,fmt) i,   layer(i)%n,&
                            SWITCH_coriolis, &
                            SWITCH_inertia
    end do
    write(200,'(A)') '  '
    write(200,'(A36,I0)')    'Number of points per period       = ', num_time_steps
    write(200,'(A36,F16.4)') 'Time step (s)                     = ', time_step
    close(200)
end subroutine write_summary


subroutine write_radial_profiles
    implicit none
    integer       :: il, ir
    character(200) :: fname, fmt

    write(fname, '(A)') trim(dir_output)//"/radial_profile_interface.dat"
    open (10, file=fname)
    fmt = "(8(A30))"
    write(10,fmt)   '#1)Radius(m)', &
                    '#2)Density(kg/m3)', &
                    '#3)Shear_modulus(Pa)', &
                    '#4)Bulk_modulus(Pa)', &
                    '#5)Viscosity(Pa s)', &
                    '#6)Density_gradient(kg/m4)', &
                    '#7)Mass_(kg)', &
                    '#8)Gravitational_acc(m/s2)'
    fmt = "(2(F30.8), 5(ES30.8), F30.8)"
    do il = 1, numoflayers
        do ir = 1, layer(il)%n
            write(10,fmt)   layer(il)%r_i(ir)    , & ! !1
                            layer(il)%rho_i(ir)  , &
                            layer(il)%mu_i(ir)   , &
                            layer(il)%K_i(ir)    , &
                            layer(il)%eta_i(ir)  , &
                            layer(il)%rho_r_i(ir), &
                            layer(il)%mass_i(ir) , &
                            layer(il)%g_i(ir)
        end do
    end do
    close(10)

    write(fname, '(A)') trim(dir_output)//"/radial_profile_layer.dat"
    open(10, file=fname)
    fmt="(5(A30))"
    write(10,fmt)   '#1)Radius(m)', &
                    '#2)Density(kg/m3)', &
                    '#3)Density_gradient(kg/m4)', &
                    '#4)Mass(kg)', &
                    '#5)Gravitational_acc(m/s2)'
    fmt="(2(F30.8), 2(ES30.8), F30.8)"
    do il = 1, numoflayers
        do ir = 1, layer(il)%n+1
            write(10,fmt)   layer(il)%r_l(ir)    , &
                            layer(il)%rho_l(ir)  , &
                            layer(il)%rho_r_l(ir), &
                            layer(il)%mass_l(ir) , &
                            layer(il)%g_l(ir)
        end do
    end do
    close(10)
end subroutine write_radial_profiles



subroutine write_spectrum(it2)
    implicit none
    integer       :: it2, j, i
    character(200) :: fname

    write(fname, '(A)') trim(dir_output)//"/spectrum.dat"
    open(400, file=fname)
    write(400,'(A, I6)') '#Time step is ', it2
    do j = 1, jmax
        write(400, '(I30, 12(ES30.16))') j, (/ (layer(i)%spectrum(j), i = 1, numoflayers), (layer(i)%spectrum_v(j), i = 1, numoflayers) /)
    end do
    close(400)
end subroutine write_spectrum

subroutine open_files_time_series
    implicit none
    integer :: i
    character(len=100) :: frmt
    character(len=500) :: header_dissipation
    iu = (/ (100+i, i = 1,20) /)

    open( iu(1), file=trim(dir_output)//"/radial_displacement_interface20_22.dat")
    write(iu(1), '(6(A20))') '# 1)Time(P)', '2)Core/Mantle', '3)Mantle/HPIce', '4)HPIce/Ocean', '5)Ocean/IceI', '6)Surface'
    
    open (iu(3), file=trim(dir_output)//"/love_numbers_time.dat")
    write(iu(3),'(7(A20), 3(A20))') ' # 1)Time(P)', '2)k20', '3)Re(k22)', '4)Im(k22)', '5)h20', '6)Re(h22)', '7)Im(h22)', '8)l20', '9)Re(l22)', '10)Im(l22)'

    open (iu(4), file=trim(dir_output)//"/dissipation.dat")
    write(iu(4),'(6(A20))') '# 1)Time(P)', '2)Core(GW)', '3)Mantle(GW)', '4)HPIce(GW)', '5)Ocean(GW)', '6)IceI(GW)'

    open (iu(5), file=trim(dir_output)//"/love_numbers_period.dat")
    write(iu(5),'(8(A20),3(A20),9(A20))') '#Time(P)', 'k20', 'Re(k22)', 'Im(k22)', 'h20', 'Re(h22)', 'Im(h22)', 'Dissipation (GW)', 'l20', 'Re(l22)', 'Im(l22)',&
    'k40', 'Re(k42)', 'Im(k42)', 'k60', 'Re(k62)', 'Im(k62)', 'k80', 'Re(k82)', 'Im(k82)'
    close(iu(5))

    open (iu(6), file=trim(dir_output)//"/peaks.dat")
    write(iu(6),'(7(A20))') '#1: k20', '2: Re(k22)', '3: Im(k22)', '4: h20', '5: Re(h22)', '6: Im(h22)', '7: Dissipation (GW)'

    open (iu(7), file=trim(dir_output)//"/time_series.dat")
    write(frmt,'(A,I0,A)') '((A40), ', numoflayers, '(I30), 3(A30))'
    write(header_dissipation, frmt) '#Dissipation (GW) layers - Time (P)', (/ (i, i = 1,numoflayers) /), 'Vg20', 'Re(Vg22)', 'Im(Vg22)'
    write(iu(7), '(A)') header_dissipation

    open(iu(8), file=trim(dir_output)//"/periodic_data.data")
    write(frmt,'(A,I0,A)') '((A40), ', numoflayers, '(I30), 4(A30))'
    write(header_dissipation, frmt) '#Dissipation (GW) layers - Period (P)', (/ (i, i = 1,numoflayers) /), 'k20', 'Re(k21)', 'Re(k22)', 'Im(k22)'
    write(iu(8), '(A)') header_dissipation
    
end subroutine open_files_time_series

subroutine open_file_full_potential
    implicit none
    iu(9) = 129
    open(iu(9), file=trim(dir_output)//"/full_potential.dat")
    write(iu(9), '(A)') '### Full gravitational potential ceofficients at the surface of the body ###'
    write(iu(9), '(A, I0, A)') '### Columns: Time(P), potential coefficinets up to degree ', 10, ' with strucuture (0, Re(2), Im(2))####'
end subroutine open_file_full_potential

subroutine write_file_full_potential(t2)
    implicit none
    real(dp) :: t2
    integer  :: jmi, il, nn, j, m
    il = numoflayers
    nn = layer(il)%n

    do j = 2, jmax_potential, 2
        do m = 0, min(j, mmax)
            jmi = jm(j, m)
            write(iu(9), '(F20.8, 2(I3), 2(ES25.12E2), 2(ES25.12E2))') t2/rotation_period, j, m, real(layer(il)%pot(jmi, nn)), imag(layer(il)%pot(jmi, nn)), &
                                                                                   real(layer(il)%ur_i(jmi, 2)), imag(layer(il)%ur_i(jmi, 2))
        end do
    end do

end subroutine write_file_full_potential


subroutine open_files_time_series_rerun
    implicit none
    integer :: i
    iu = (/ (100+i, i = 1,20) /)

    open( iu(1), file=trim(dir_output)//"/radial_displacement_interface20_22_r.dat")
    write(iu(1), '(6(A20))') '# 1)Time(P)', '2)Core/Mantle', '3)Mantle/HPIce', '4)HPIce/Ocean', '5)Ocean/IceI', '6)Surface'
    
    open (iu(3), file=trim(dir_output)//"/love_numbers_time_r.dat")
    write(iu(3),'(7(A20), 3(A20))') ' # 1)Time(P)', '2)k20', '3)Re(k22)', '4)Im(k22)', '5)h20', '6)Re(h22)', '7)Im(h22)', '8)l20', '9)Re(l22)', '10)Im(l22)'

    open (iu(4), file=trim(dir_output)//"/dissipation_r.dat")
    write(iu(4),'(6(A20))') '# 1)Time(P)', '2)Core(GW)', '3)Mantle(GW)', '4)HPIce(GW)', '5)Ocean(GW)', '6)IceI(GW)'

    ! open (iu(5), file=trim(dir_output)//"/love_numbers_period.dat")
    ! write(iu(5),'(8(A20),3(A20),9(A20))') '#Time(P)', 'k20', 'Re(k22)', 'Im(k22)', 'h20', 'Re(h22)', 'Im(h22)', 'Dissipation (TW)', 'l20', 'Re(l22)', 'Im(l22)',&
    ! 'k40', 'Re(k42)', 'Im(k42)', 'k60', 'Re(k62)', 'Im(k62)', 'k80', 'Re(k82)', 'Im(k82)'
    ! close(iu(5))
end subroutine open_files_time_series_rerun

subroutine write_time_series(it2, t2)
    implicit none
    integer  :: it2, il, nn
    real(dp) :: t2
    character(len=200) :: frmt
    !! The dissipation as a function of time
    il = numoflayers; nn =layer(il)%n
    write(frmt,'(A,I0,A)') '((F40.12), ', numoflayers, '(F30.12), 3(F30.12))'
    write(iu(7), frmt) t2/rotation_period, layer(:)%dissipation*W_to_GW,    real(layer(il)%pot(jm(2,0), nn)), &
                                                                            real(layer(il)%pot(jm(2,2), nn)), &
                                                                            imag(layer(il)%pot(jm(2,2), nn))

    if(mod(it2, num_time_steps) == 0) then
        write(frmt,'(A,I0,A)') '((F40.12), ', numoflayers, '(F30.12), 4(F30.6))'
        write(iu(8), frmt) t2/rotation_period, layer(:)%average_dissip*W_to_GW, real(love_k(jm(2,0))), &
                                                                                real(love_k(jm(2,1))), &
                                                                                real(love_k(jm(2,2))), &
                                                                                imag(love_k(jm(2,2)))
        layer(:)%average_dissip = 0.5_dp * layer(:)%dissipation
    end if
end subroutine write_time_series



subroutine write_time_series_rerun(t2)
    implicit none
    real(dp) :: t2
    !!j = 2 m = 0
    ! write(iu(1), '(4(F20.8))') t2/rotation_period,  real( layer(1)%ur_i(layer(1)%n, jm2(2,0))), &
    !                                                 real( layer(2)%ur_i(layer(2)%n, jm2(2,0))), &
    !                                                 real( layer(3)%ur_i(layer(3)%n, jm2(2,0)))

    ! !!j = 2 m = 2 - real
    ! write(iu(1), '(4(F20.8))') t2/rotation_period,  real( layer(1)%ur_i(layer(1)%n, jm2(2,2))), &
    !                                                 real( layer(2)%ur_i(layer(2)%n, jm2(2,2))), &
    !                                                 real( layer(3)%ur_i(layer(3)%n, jm2(2,2)))
    ! !!j = 2 m = 0 - imag
    ! write(iu(1), '(4(F20.8))') t2/rotation_period,  aimag(layer(1)%ur_i(layer(1)%n, jm2(2,2))), &
    !                                                 aimag(layer(2)%ur_i(layer(2)%n, jm2(2,2))), &
    !                                                 aimag(layer(3)%ur_i(layer(3)%n, jm2(2,2)))
    ! !! Gravitational potential as a function of time
    ! write(iu(3), '(10(F20.8))') t2/rotation_period, k20, real_k22, imag_k22, &
    !                                                 h20, real_h22, imag_h22, &
    !                                                 l20, real_l22, imag_l22
    !! The dissipation as a function of time
    write(iu(4), '(4(F20.8))') t2/rotation_period, layer(:)%dissipation*W_to_GW
end subroutine write_time_series_rerun

subroutine close_files
    implicit none
    
    close(iu(1))
    close(iu(2))
    close(iu(3))
    close(iu(4))
    close(iu(7))
    close(iu(8))
end subroutine close_files

subroutine write_tidal_potential
    implicit none
    integer :: il, nn, it
    real(dp) :: tt
    il = numoflayers; nn = layer(il)%n
    open(10, file=trim(dir_output)//"/tidal_potential.dat")
    if(SWITCH_mmtides == 1) write(10,'(8(A20))') '#Time(s)', 'Time(P)', 'Ecc_V20', 'Ecc_RV22', 'Ecc_IV22', 'MM_V20', 'MM_RV22', 'MM_IV22'
    if(SWITCH_mmtides == 0) write(10,'(5(A20))') '#Time(s)', 'Time(P)', 'Ecc_V20', 'Ecc_RV22', 'Ecc_IV22'

    do it = 1, num_time_steps+1
        tt = (it-1)*time_step
        ! if(SWITCH_mmtides == 1) then
        ! write(10,'(8(F20.8))') tt, tt/rotation_period,  real(tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 0, tt, 'p')), &
        !                                                 real(tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p')), &
        !                                                 imag(tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p')), &
        !                                                 real(mm_tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 0, tt, 'p')), &
        !                                                 real(mm_tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p')), &
        !                                                 imag(mm_tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p'))
        ! elseif(SWITCH_mmtides == 2) then
        !     write(10,'(5(F20.8))') tt, tt/rotation_period,  real(lunar_tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 0, tt, 'p')), &
        !                                                     real(lunar_tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p')), &
        !                                                     imag(lunar_tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p'))
        ! else
        !     write(10,'(5(F20.8))') tt, tt/rotation_period,  real(tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 0, tt, 'p')), &
        !                                                     real(tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p')), &
        !                                                     imag(tides(layer(il)%rho_i(nn), layer(il)%r_i(nn), 2, tt, 'p'))
        ! end if
    end do
    close(10)
end subroutine write_tidal_potential




subroutine write_initial_signal
    implicit none
    signal = 1
    
    open(400, file=trim(dir_output)//"/signal")
    write(400, *) signal
    close(400)
    
end subroutine write_initial_signal


subroutine write_density_f(t2, it2)
    implicit none
    real(dp) :: t2
    integer :: it2, il
    character(300) :: fname
    
    write(fname, '(A, I0)') trim(dir_output)//'/sh_coeffs/density/density.', (it2 - start_output*num_time_steps)/write_output
    open(10, file=fname, form='unformatted',access='stream',status='replace')
    write(10) jmax, t2, numoflayers
    write(10) layer(:)%n
    do il = 1, numoflayers
        write(10) layer(il)%rhoe
    end do
    close(10)
end subroutine write_density_f

subroutine write_velocity_f(t2, it2)
    implicit none
    real(dp) :: t2
    integer :: it2, il
    character(300) :: fname
    
    write(fname, '(A, I0)') trim(dir_output)//'/sh_coeffs/velocity/velocity.', (it2 - start_output*num_time_steps)/write_output
    open(10, file=fname, form='unformatted',access='stream',status='replace')
    write(10) jmax, t2, numoflayers
    write(10) layer(:)%n
    do il = 1, numoflayers
        write(10) layer(il)%vel
    end do
    close(10)
end subroutine write_velocity_f



subroutine write_stress_f(t2, it2)
    implicit none
    real(dp) :: t2
    integer :: it2, il
    character(300) :: fname

    write(fname, '(A, I0)') trim(dir_output)//'/sh_coeffs/stress/stress.', (it2 - start_output*num_time_steps)/write_output
    open(10, file=fname, form='unformatted',access='stream',status='replace')
    write(10) jmax, t2, numoflayers
    write(10) layer(:)%n
    do il = 1, numoflayers
        write(10) layer(il)%str
    end do
    close(10)
end subroutine write_stress_f

subroutine write_displacement_f(t2, it2)
    implicit none
    real(dp) :: t2
    integer :: it2, il
    character(300) :: fname

    write(fname, '(A, I0)') trim(dir_output)//'/sh_coeffs/displacement/displacement.', (it2 - start_output*num_time_steps)/write_output
    open(10, file=fname, form='unformatted',access='stream',status='replace')
    write(10) jmax, t2, numoflayers
    write(10) layer(:)%n
    do il = 1, numoflayers
        write(10) layer(il)%dis
    end do
    close(10)
end subroutine write_displacement_f

subroutine write_checkpoint(t2, it2)
    implicit none
    real(dp) :: t2
    integer :: it2, il
    character(300) :: fname, fname_zip
    write(fname, '(A, I0)') trim(dir_output)//'/checkpoint/checkpoint.', it2
    open(10, file=fname,form='unformatted',access='stream',status='replace')
    write(10) jmax, t2
    do il = 1, numoflayers
        write(10) layer(il)%rhoe
    end do
    do il = 1, numoflayers
        write(10) layer(il)%vel
    end do
    do il = 1, numoflayers
        write(10) layer(il)%pvel
    end do
    do il = 1, numoflayers
        write(10) layer(il)%str
    end do
    do il = 1, numoflayers
        write(10) layer(il)%dis
    end do
    close(10)
    write(fname_zip, '(A)') trim(dir_output)//'/checkpoint/cp.zip'
    if(it2 == checkpoint)then
        call system("zip -j "//trim(fname_zip)//" "//trim(fname))
    else
        call system("zip -uj "//trim(fname_zip)//" "//trim(fname))
    end if
    call system("rm -f "//trim(fname))
    write(fname, '(A)') trim(dir_output)//'/checkpoint/last_checkpoint_number'
    open(10, file=fname)
    write(10, *) it2
    close(10)
end subroutine write_checkpoint

end module mod_write_f