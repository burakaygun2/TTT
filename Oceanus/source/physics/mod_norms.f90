module mod_norms
    use mod_parameters
    use mod_variables
    use mod_write_f
    use mod_write_h5
    use mod_math
    implicit none
    
contains
    

subroutine compute_dissipation(it2)
    implicit none
    integer :: i, k, j, m, jmi
    real(dp) :: dr
    real(dp) :: r2_0, r2_1, e0, e1, str0, str1, lay_int, val
    integer  :: it2
    layer(:)%dissipation = 0.0d0
    do i = 1, numoflayers
        layer(i)%spectrum = 0.0d0
    end do
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP PRIVATE(jmi, j, m, i, k, dr, r2_0, r2_1, e0, e1, str0, str1, lay_int, val)
    !$OMP DO SCHEDULE(dynamic,1)
    do jmi = 2, compute_scalar_harm
        j = jmindx(jmi, 2); m = jmindx(jmi, 3)
        do i = 1, numoflayers
            lay_int = 0.0_dp
            do k = 1, layer(i)%n - 1
                dr   = layer(i)%r_i(k+1) - layer(i)%r_i(k)
                r2_0 = layer(i)%r_i(k  )**2
                r2_1 = layer(i)%r_i(k+1)**2
                e0   = 1.0_dp / (2.0_dp*layer(i)%eta_i(k  )*hyperviscosity(j))
                e1   = 1.0_dp / (2.0_dp*layer(i)%eta_i(k+1)*hyperviscosity(j))
                str0   =    abs(layer(i)%str(6*(jmi-1)+2, k))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+3, k))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+4, k))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+5, k))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+6, k))**2
                            
                str1   =    abs(layer(i)%str(6*(jmi-1)+2, k+1))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+3, k+1))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+4, k+1))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+5, k+1))**2 + &
                            abs(layer(i)%str(6*(jmi-1)+6, k+1))**2
                lay_int = lay_int + 0.5_dp * dr * (str0*r2_0*e0 + str1*r2_1*e1)
            end do
            val = lay_int * real(2 - dirac(m, 0), dp)
            !$OMP ATOMIC
            layer(i)%dissipation = layer(i)%dissipation   + val
            !$OMP ATOMIC
            layer(i)%spectrum(j) = layer(i)%spectrum(j)   + val
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    do i = 1, numoflayers
        layer(i)%dissipation_time(mod(it2-1, num_time_steps) + 1) = layer(i)%dissipation
    end do
    layer(:)%average_dissip = layer(:)%average_dissip + layer(:)%dissipation
    if(mod(it2, num_time_steps) == 0) then
        layer(:)%average_dissip = (layer(:)%average_dissip - layer(:)%dissipation*0.5_dp)/real(num_time_steps, dp)
    end if
    if (mod(it2, 20) == 0) call check_spectrum(it2)
end subroutine compute_dissipation


subroutine velocity_norm
    implicit none
    integer    :: il, ir, j, m, nn, jmi
    real(dp)    :: integral(1:numoflayers)
    real(dp)    :: vv0, r2_0, r2_1
    real(dp)    :: vv1, dr

    do il = 1 , numoflayers
        layer(il)%v_norm = 0.0d0
        layer(il)%spectrum_v = zero
    end do
    do jmi = 2, compute_vector_harm, 3
        j = jmlindx(jmi, 2); m = jmlindx(jmi, 3); 
        do il = 1, numoflayers
            integral(il) = zero
            nn = layer(il)%n
            do ir = 2, nn
                dr = layer(il)%r_l(ir+1) - layer(il)%r_l(ir)
                r2_0 = layer(il)%r_l(ir  )**2
                r2_1 = layer(il)%r_l(ir+1)**2
                vv0 =   abs(layer(il)%vel(jmi  , ir))**2 + &
                        abs(layer(il)%vel(jmi+1, ir))**2 + &
                        abs(layer(il)%vel(jmi+2, ir))**2

                vv1 =   abs(layer(il)%vel(jmi  , ir+1))**2 + &
                        abs(layer(il)%vel(jmi+1, ir+1))**2 + &
                        abs(layer(il)%vel(jmi+2, ir+1))**2

                integral(il) = integral(il) + dr * 0.5_dp * (vv0*r2_0 + vv1*r2_1)
            end do
            integral(il)=integral(il)* real(2-dirac(m,0),dp)
            layer(il)%v_norm = layer(il)%v_norm   + integral(il)
            layer(il)%spectrum_v(j) = layer(il)%spectrum_v(j) + integral(il)
        end do
    end do
    ! do i = 1, numoflayers
    !     layer(i)%dissipation_time(mod(it2-1, num_time_steps) + 1) = layer(i)%dissipation
    ! end do
    ! layer(:)%average_dissip = layer(:)%average_dissip + layer(:)%dissipation
    ! if(mod(it2, num_time_steps) == 0) then
    !     layer(:)%average_dissip = (layer(:)%average_dissip - layer(:)%dissipation*0.5_dp)/real(num_time_steps, dp)
    ! end if
end subroutine velocity_norm


subroutine check_spectrum(it2)
    implicit none
    integer :: i, change_status, it2
    real(dp) :: trunc_3

    ! call write_spectrum(it2)
    if(mod(it2, num_time_steps)==0)  call write_spectrum_h5
    change_status = 0
    old_jmax = jmax
    if((jmax == Max_jmax).or.(jmax==2)) goto 8
    do i = 1, numoflayers
        if(maxval(layer(i)%spectrum_v(3:6)) > 0.0_dp) trunc_3 = log10(maxval(layer(i)%spectrum_v(3:6))) - 4.0_dp
        if((jmax.ge.3) .and. minval(layer(i)%spectrum_v(jmax-5:jmax))>0.0_dp) then
            if(maxval(log10(layer(i)%spectrum_v(jmax-5:jmax))).gt.trunc_3) then
                change_status = 1
                exit
            end if
        end if
    end do
    if(change_status.eq.1)then
        jmax = jmax + 4
        if(jmax > Max_jmax) jmax = Max_jmax
        jfinal = jmax
        if (SWITCH_vgradv == 1) then
            compute_scalar_harm = jm (jmax, jmax)
            compute_vector_harm = jml(jmax, jmax, jmax+1)
        else
            compute_scalar_harm = jm (jmax, mmax)
            compute_vector_harm = jml(jmax, mmax, jmax+1)
        end if
    end if
8 continue
end subroutine check_spectrum

subroutine overwritejmax(lmax)
    implicit none
    integer :: lmax
    jmax = lmax
    jfinal = jmax
    if (SWITCH_vgradv == 1) then
        compute_scalar_harm = jm (jmax, mmax)
        compute_vector_harm = jml(jmax, mmax, jmax+1)
    else
        compute_scalar_harm = jm (jmax, mmax)
        compute_vector_harm = jml(jmax, mmax, jmax+1)
    end if
end subroutine overwritejmax

subroutine check_stability(p)
    implicit none
    integer :: p, i
    real(dp) :: delta_k22   (1:9)
    real(dp) :: diff1, diff2
    do i = 1, 10-1
        delta_k22(i) = abs(rk22_p(p-i+1) -rk22_p(p-i))
        delta_k22(i) = delta_k22(i)*100.0_dp/rk22_p(p-i)
    end do
    if (SWITCH_coriolis == 0) then
        diff1 = abs(rk22_p(p) - k20_p(p))
        diff2 = abs(rk22_p(p) - ik22_p(p))
        if ((maxval(delta_k22) <= 0.01d0).and.(diff1<1d-6).and.(diff2<1d-6)) signal=0
    else
        if (maxval(delta_k22) <= 0.01d0) signal=0
    end if

end subroutine check_stability

subroutine check_signal
    implicit none
    
    if(signal==1)then
        open(400, file=trim(dir_output)//"/signal")
        read(400, *) signal
        close(400)
    end if
end subroutine check_signal

end module mod_norms