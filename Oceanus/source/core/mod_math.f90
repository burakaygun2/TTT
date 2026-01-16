module mod_math
    use mod_parameters
    use mod_variables
    use mod_matrix_coeffs
    use mod_potential
    implicit none
    type :: div_s_c
        real(dp) :: jm1_dsj0, jm1_sj0, jm1_dsj2, jm1_sj2, jm1_dsjm2, jm1_sjm2
        real(dp) :: jp1_dsj0, jp1_sj0, jp1_dsj2, jp1_sj2, jp1_dsjp2, jp1_sjp2
        real(dp) :: jm0_dsjm1, jm0_sjm1, jm0_dsjp1, jm0_sjp1
        real(dp) :: jm1_rho, jp1_rho
    end type
    type :: grad_v_c
        real(dp) :: jm2_dvjm1, jm2_vjm1, j2_dvjm1, j2_vjm1, j2_dvjp1, j2_vjp1, jp2_dvjp1, jp2_vjp1
    end type
    type(div_s_c)  :: dsc
    type(grad_v_c) :: gvc
contains

subroutine peaks(tt, it2)
    implicit none
    real(dp) :: tt
    integer :: it2
    real(dp) :: Delta_20_1 , Delta_20_2
    real(dp) :: Delta_r22_1, Delta_r22_2
    real(dp) :: Delta_i22_1, Delta_i22_2
    real(dp) :: Delta_d_1  , Delta_d_2
    real(dp) :: fp
    real(dp) :: diff1, diff2


    if(it2 == 3) consecutive_count = 0
    if(it2<3) go to 10
    gp20_2 = k20; rgp22_2 = real_k22; igp22_2 = imag_k22;
    ! dissipation_2 = layer(2)%dissipation*W_to_GW

    Delta_20_1 = gp20_2 - gp20_1
    Delta_20_2 = gp20_1 - gp20_0
    if((Delta_20_1<0).and.(Delta_20_2>0)) then
        fp =  ang_vel * ang_vel * ecc * sqrt(18.0_dp * pi/10.0_dp) * (rad(numoflayers)**2)
        write(iu(6), '(I5, F20.5, F20.5)') 1, (tt - time_step)/rotation_period, gp20_1/fp
    end if

    Delta_r22_1 = rgp22_2 - rgp22_1
    Delta_r22_2 = rgp22_1 - rgp22_0
    if((Delta_r22_1<0).and.(Delta_r22_2>0)) then
        peaks_r22(3) = rgp22_1
        peak_time_r22(3) = tt - time_step
        fp =  ang_vel * ang_vel * ecc * sqrt(27.0_dp * pi/10.0_dp) * (rad(numoflayers)**2)

        write(iu(6), '(I5, F20.5, F20.5)') 2, peak_time_r22(3)/rotation_period, peaks_r22(3)/fp
        if((abs(peaks_r22(3)-peaks_r22(2)) < 1d-5).and.((abs(abs(peak_time_r22(3)-peak_time_r22(2)) - rotation_period)))<1d-1) then
            if (SWITCH_coriolis == 1) then
                consecutive_count = consecutive_count + 1
                if(consecutive_count >= 5) then
                    signal = 0
                    print*, 'Peak to peak stability has been reached...'
                end if
            else
                consecutive_count = consecutive_count + 1
                diff1 = abs(igp22_1/fp - peaks_r22(3)/fp)
                diff2 = abs(gp20_1/fp  - peaks_r22(3)/fp)
                if ((consecutive_count >= 5).and.(diff1 < 1d-5).and.(diff2 < 1d-5)) then
                    signal = 0
                    print*, 'Peak to peak stability has been reached...'
                end if
            end if
        end if
        peaks_r22(1) = peaks_r22(2)
        peaks_r22(2) = peaks_r22(3)
        peak_time_r22(1) = peak_time_r22(2)
        peak_time_r22(2) = peak_time_r22(3)
    end if
    
    Delta_i22_1 = igp22_2 - igp22_1
    Delta_i22_2 = igp22_1 - igp22_0
    if((Delta_i22_1<0).and.(Delta_i22_2>0)) then
        fp = ang_vel * ang_vel * ecc * sqrt(48.0_dp * pi/10.0_dp) * (rad(numoflayers)**2)
        write(iu(6), '(I5, F20.5, F20.5)') 3, (tt - time_step)/rotation_period, igp22_1/fp
    end if
    
    Delta_d_1 = dissipation_2 - dissipation_1
    Delta_d_2 = dissipation_1 - dissipation_0
    
    if((Delta_d_1<0).and.(Delta_d_2>0)) then
        ! print*, 'Dissipaiton peak at...', (tt-time_step)/rotation_period, dissipation_1
        pt_d_1 = pt_d_2
        pt_d_2 = tt
        ! print*, 'Peak time difference (dissipation)...', (pt_d_2-pt_d_1)/rotation_period
    end if

    gp20_0 = gp20_1
    gp20_1 = gp20_2


    rgp22_0 = rgp22_1
    rgp22_1 = rgp22_2


    igp22_0 = igp22_1
    igp22_1 = igp22_2

    dissipation_0 = dissipation_1
    dissipation_1 = dissipation_2
10 continue
end subroutine peaks


subroutine determine_dt
    implicit none
    integer :: j, m, ir, il, oceindx
    integer :: jm1, jp1, jmind
    real(dp) :: vr, vt, erjm1, erjp1, jj, etjm1, etjp1
    real(dp) :: cr, ct, dr, cd

    do jmind = 2, compute_scalar_harm
        j = jmindx(jmind, 2); m = jmindx(jmind, 3); jj = real(j, dp);
        jm1 = jml(j, m, j-1); jp1=jml(j, m, j+1)
        erjm1 =  sqrt( jj       /(2.0_dp*jj+1.0_dp))
        erjp1 = -sqrt((jj+1.0_dp)/(2.0_dp*jj+1.0_dp))
        etjm1 =  (jj+1.0_dp)           / (2.0_dp*jj+1.0_dp)
        etjp1 =  sqrt(jj*(jj+1.0_dp)) / (2.0_dp*jj+1.0_dp)
        do il = 1, num_oceans
            oceindx = layer_ocean(il)
            do ir = 2, layer(il)%n
                dr = abs(layer(oceindx)%r_i(ir)-layer(oceindx)%r_i(ir-1))
                vr = abs(erjm1 * layer(oceindx)%vel(jm1, ir) + erjp1 * layer(oceindx)%vel(jp1, ir))
                vt = abs(etjm1 * layer(oceindx)%vel(jm1, ir) + etjp1 * layer(oceindx)%vel(jp1, ir))

                ! print*, time_step, vr, dr
                cd = time_step * (layer(oceindx)%eta_i(1)/layer(oceindx)%rho_i(1))/ dr**2
                cr =  time_step * (vr / dr) !**(-4.0/3.0)
                ct = vt * time_step / (layer(oceindx)%r_l(ir)/sqrt(real(jmax)*(real(jmax)+1)))
                if ((cr>1.0_dp) .or. (ct>1.0_dp) .or. (cd>1.0_dp)) then
                    print*, j, m, cd, cr, ct, vr, vt
                    stop
                    if(dr / vr < dtmin_r) then
                        dtmin_r = (dr / vr) ** (4.0 / 3.0)
                        change_time_step = 0
                    end if
                end if
            end do
        end do
    end do
end subroutine determine_dt


subroutine Divergence_tensor
    implicit none
    integer :: j, m, il, ir, nn
    integer :: jm1, jm0, jp1, jms
    complex(dp) :: dsj0, asj0, dsj2, asj2, dsjm2, asjm2, dsjp2, asjp2, rho_jm1, rho_jp1, dsjm1, asjm1, dsjp1, asjp1
    real(dp)    :: hi, hup2, hlo2, hup2_r, hlo2_r
    
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn-1
            hi   =   layer(il)%r_i(ir+1)   - layer(il)%r_i(ir  )
            hup2 =  (layer(il)%r_i(ir+1)   - layer(il)%r_l(ir+1)) / hi;           hup2_r = hup2/layer(il)%r_l(ir+1)
            hlo2 =  (layer(il)%r_l(ir+1)   - layer(il)%r_i(ir  )) / hi;           hlo2_r = hlo2/layer(il)%r_l(ir+1)
            do j = jstart, jfinal
                call div_tensor_coeffs(real(j, dp))
                do m = 0, min(2,j), 2
                    jm1 = jml2(j,m,j-1)
                    jm0 = jml2(j,m,j  )
                    jp1 = jml2(j,m,j+1)
                    jms = (jm2(j,m) - 1)*6

                    dsj0  = dsc%jm1_dsj0  * (layer(il)%str(ir+1, jms+1)        - layer(il)%str(ir, jms+1)) / hi
                    asj0  = dsc%jm1_sj0   * (layer(il)%str(ir+1, jms+1)*hlo2_r + layer(il)%str(ir, jms+1)*hup2_r)
                    dsjm2 = dsc%jm1_dsjm2 * (layer(il)%str(ir+1, jms+2)        - layer(il)%str(ir, jms+2)) / hi
                    asjm2 = dsc%jm1_sjm2  * (layer(il)%str(ir+1, jms+2)*hlo2_r + layer(il)%str(ir, jms+2)*hup2_r)
                    dsj2  = dsc%jm1_dsj2  * (layer(il)%str(ir+1, jms+4)        - layer(il)%str(ir, jms+4)) / hi
                    asj2  = dsc%jm1_sj2   * (layer(il)%str(ir+1, jms+4)*hlo2_r + layer(il)%str(ir, jms+4)*hup2_r)
                    rho_jm1 = dsc%jm1_rho * (layer(il)%g_l(ir+1) * (layer(il)%rhoe(ir+1, jm2(j,m))*hlo2 + layer(il)%rhoe(ir, jm2(j,m))*hup2))

                    layer(il)%div_tensor(ir, jm1) =  dsj0 +  asj0 + dsj2 + asj2 + dsjm2 + asjm2 + rho_jm1

                    dsj0  = dsc%jp1_dsj0  * (layer(il)%str(ir+1, jms+1)        - layer(il)%str(ir, jms+1)) / hi
                    asj0  = dsc%jp1_sj0   * (layer(il)%str(ir+1, jms+1)*hlo2_r + layer(il)%str(ir, jms+1)*hup2_r)
                    dsjp2 = dsc%jp1_dsjp2 * (layer(il)%str(ir+1, jms+6)        - layer(il)%str(ir, jms+6)) / hi
                    asjp2 = dsc%jp1_sjp2  * (layer(il)%str(ir+1, jms+6)*hlo2_r + layer(il)%str(ir, jms+6)*hup2_r)
                    dsj2  = dsc%jp1_dsj2  * (layer(il)%str(ir+1, jms+4)        - layer(il)%str(ir, jms+4)) / hi
                    asj2  = dsc%jp1_sj2   * (layer(il)%str(ir+1, jms+4)*hlo2_r + layer(il)%str(ir, jms+4)*hup2_r)
                    rho_jp1 = dsc%jp1_rho * (layer(il)%g_l(ir+1) * (layer(il)%rhoe(ir+1, jm2(j,m))*hlo2 + layer(il)%rhoe(ir, jm2(j,m))*hup2))

                    layer(il)%div_tensor(ir, jp1) =  dsj0 +  asj0 + dsj2 + asj2 + dsjp2 + asjp2 + rho_jp1

                    dsjm1 = dsc%jm0_dsjm1 * (layer(il)%str(ir+1, jms+3)        - layer(il)%str(ir, jms+3)) / hi
                    asjm1 = dsc%jm0_sjm1  * (layer(il)%str(ir+1, jms+3)*hlo2_r + layer(il)%str(ir, jms+3)*hup2_r)
                    dsjp1 = dsc%jm0_dsjp1 * (layer(il)%str(ir+1, jms+5)        - layer(il)%str(ir, jms+5)) / hi
                    asjp1 = dsc%jm0_sjp1  * (layer(il)%str(ir+1, jms+5)*hlo2_r + layer(il)%str(ir, jms+5)*hup2_r)

                    layer(il)%div_tensor(ir, jm0) =  dsjm1 + asjm1 + dsjp1 + asjp1
                end do
            end do
        end do
    end do
end subroutine Divergence_tensor

subroutine Grad_velocity
    implicit none
    integer    :: j, m, il, ir, nn
    integer    :: jm1, jm0, jp1, jms
    real(dp)    :: hl, hup, hlo, hup_r, hlo_r
    complex(dp) :: dvjm1, avjm1, dvjp1, avjp1


    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn
            hl  =  layer(il)%r_l(ir+1) - layer(il)%r_l(ir)
            hup = (layer(il)%r_l(ir+1) - layer(il)%r_i(ir)) / hl;           hup_r = hup/layer(il)%r_i(ir)
            hlo = (layer(il)%r_i(ir)   - layer(il)%r_l(ir)) / hl;           hlo_r = hlo/layer(il)%r_i(ir)
            do j = jstart, jfinal
                call grad_v_coeffs(real(j, dp))
                do m = 0, min(2, j), 2
                    jm1 = jml2(j,m,j-1)
                    jm0 = jml2(j,m,j  )
                    jp1 = jml2(j,m,j+1)
                    jms = (jm2(j,m) - 1)*6

                    dvjm1 = gvc%jm2_dvjm1 * (layer(il)%vel(ir+1, jm1)       - layer(il)%vel(ir, jm1)) / hl
                    avjm1 = gvc%jm2_vjm1  * (layer(il)%vel(ir+1, jm1)*hlo_r + layer(il)%vel(ir, jm1)*hup_r)
                    layer(il)%grad_v(ir, jms+2) = dvjm1 + avjm1

                    dvjm1 = gvc%j2_dvjm1 * (layer(il)%vel(ir+1, jm1) - layer(il)%vel(ir, jm1)) / hl
                    avjm1 = gvc%j2_vjm1  * (layer(il)%vel(ir+1, jm1)*hlo_r + layer(il)%vel(ir, jm1)*hup_r)
                    dvjp1 = gvc%j2_dvjp1 * (layer(il)%vel(ir+1, jp1) - layer(il)%vel(ir, jp1)) / hl
                    avjp1 = gvc%j2_vjp1  * (layer(il)%vel(ir+1, jp1)*hlo_r + layer(il)%vel(ir, jp1)*hup_r)
                    layer(il)%grad_v(ir, jms+4) = dvjm1 + avjm1 + dvjp1 + avjp1

                    dvjm1 = gvc%jp2_dvjp1 * (layer(il)%vel(ir+1, jp1) - layer(il)%vel(ir, jp1)) / hl
                    avjp1 = gvc%jp2_vjp1  * (layer(il)%vel(ir+1, jp1)*hlo_r + layer(il)%vel(ir, jp1)*hup_r)
                    layer(il)%grad_v(ir, jms+6) = dvjm1 + avjp1

                end do
            end do
        end do
    end do
end subroutine Grad_velocity


subroutine div_tensor_coeffs(j)
    implicit none
    real(dp) :: j

    ! !Y_{JM}^{J-1}
    dsc%jm1_dsj0  = -sqrt( j / (3.0_dp * (2.0_dp * j + 1.0_dp)))
    dsc%jm1_dsjm2 =  sqrt((j - 1.0_dp) / (2.0_dp * j - 1.0_dp))
    dsc%jm1_dsj2  = -sqrt( (j + 1.0_dp)*(2.0_dp * j + 3.0_dp) /(6._dp*(2._dp*j - 1._dp)*(2._dp*j + 1._dp) ))
    dsc%jm1_sj0   =  dsc%jm1_dsj0  * (j + 1.0_dp)
    dsc%jm1_sjm2  = -dsc%jm1_dsjm2 * (j - 2.0_dp)
    dsc%jm1_sj2   =  dsc%jm1_dsj2  * (j + 1.0_dp)
    dsc%jm1_rho   = -sqrt( j       /(2.0_dp * j + 1.0_dp))

    ! !Y_{JM}^{J+1}
    dsc%jp1_dsj0   =  sqrt((j + 1.0_dp)/((3._dp*(2._dp*j + 1._dp))))
    dsc%jp1_dsjp2  = -sqrt((j +2._dp) / (2._dp*j + 3._dp))
    dsc%jp1_dsj2   =  sqrt(( j*(2._dp*j - 1._dp) )/(6._dp*(2._dp*j + 3._dp)*(2._dp*j + 1._dp)))
    dsc%jp1_sj0    = -dsc%jp1_dsj0  * j
    dsc%jp1_sjp2   =  dsc%jp1_dsjp2 * (j + 3.0_dp)
    dsc%jp1_sj2    = -dsc%jp1_dsj2  * j
    dsc%jp1_rho    = -(-sqrt((j+1.0_dp)/(2.0_dp * j + 1.0_dp)))


    dsc%jm0_dsjm1 =  sqrt((j - 1.0_dp)/(2.0_dp*(2.0_dp * j + 1.0_dp)))
    dsc%jm0_sjm1  = -(j - 1.0_dp) * dsc%jm0_dsjm1
    dsc%jm0_dsjp1 = -sqrt((j + 2.0_dp)/(2.0_dp*(2.0_dp * j + 1.0_dp)))
    dsc%jm0_sjp1  =  (j + 2.0_dp) * dsc%jm0_dsjp1

end subroutine div_tensor_coeffs

subroutine grad_v_coeffs(j)
    implicit none
    real(dp) :: j

    gvc%jm2_dvjm1 = -sqrt((j - 1.0_dp) / (2.0_dp * j - 1.0_dp))
    gvc%jm2_vjm1  =  gvc%jm2_dvjm1 * j

    gvc%j2_dvjm1 =  sqrt( ((j + 1._dp)*(2._dp*j + 3._dp))/(6._dp*(2._dp*j - 1._dp)*(2._dp*j + 1._dp)))
    gvc%j2_dvjp1 = -sqrt((j*(2._dp*j - 1._dp)) / (6._dp*(2._dp*j + 3._dp)*(2._dp*j + 1._dp)))
    gvc%j2_vjm1  = -gvc%j2_dvjm1 * (j - 1.0_dp)
    gvc%j2_vjp1  =  gvc%j2_dvjp1 * (j + 2.0_dp)

    gvc%jp2_dvjp1 =  sqrt( (j + 2.0_dp)  / (2.0_dp * j + 3.0_dp) )
    gvc%jp2_vjp1  = -(j + 1.0_dp) * gvc%jp2_dvjp1
end subroutine grad_v_coeffs


function QST_to_vjml(j) result(Q2v)
    implicit none
    real(dp), intent(in) :: j
    real(dp), dimension(1:3) :: Q2v

    Q2v(1) =  sqrt( j        * (2.0 * j + 1.0)) / (2.0 * j + 1.0)
    Q2v(2) = -sqrt((j + 1.0) * (2.0 * j + 1.0)) / (2.0 * j + 1.0)
    Q2v(3) =  sqrt((j      ) * (      j + 1.0))
end function QST_to_vjml

subroutine find_boundary_layer_thickness(it2)
    implicit none
    integer                  :: jmind, i, il, ir, it2
    complex(dp), allocatable :: dv(:,:), dvmax(:), dl(:)
    real(dp)                 :: dr

    do i = 1, num_oceans
        il = layer_ocean(i)
        allocate(dv(compute_vector_harm, layer(il)%n))
        do ir = 1, layer(il)%n
            do jmind = 2, compute_vector_harm, 3
                dr = layer(il)%r_l(ir+1) - layer(il)%r_l(ir)
                dv(jmind  , ir) = (layer(il)%vel(jmind  , ir+1) - layer(il)%vel(jmind  , ir)) / dr
                dv(jmind+1, ir) = (layer(il)%vel(jmind+1, ir+1) - layer(il)%vel(jmind+1, ir)) / dr
                dv(jmind+2, ir) = (layer(il)%vel(jmind+2, ir+1) - layer(il)%vel(jmind+2, ir)) / dr
            end do
            ! deallocate(dv)
        end do
    end do
    allocate(dvmax(compute_vector_harm), dl(compute_vector_harm))
    do jmind = 2, compute_vector_harm!, 3
        ! dvmax(jmind  ) = cmplx(maxval(abs(real(dv(jmind  , 1:50)))), maxval(abs(imag(dv(jmind  , 1:50)))))
        ! dvmax(jmind+1) = cmplx(maxval(abs(real(dv(jmind+1, 1:50)))), maxval(abs(imag(dv(jmind+1, 1:50)))))
        ! dvmax(jmind+2) = cmplx(maxval(abs(real(dv(jmind+2, 1:50)))), maxval(abs(imag(dv(jmind+2, 1:50)))))
    end do
    do i = 1, num_oceans
        il = layer_ocean(i)
        do jmind = 2, compute_vector_harm
            do ir = 1, layer(il)%n-50
                if ( abs(real(dv(jmind  , ir)))<real(dvmax(jmind))*0.05 ) then
                    dl(jmind) = (layer(il)%r_i(ir) - layer(il)%r_i(1))
                    goto 1                
                end if
1               if  ( abs(imag(dv(jmind  , ir)))<imag(dvmax(jmind))*0.05 ) then
                    dl(jmind) = dl(jmind) + (layer(il)%r_i(ir) - layer(il)%r_i(1))*Im
                    goto 2
                end if
            end do
2           continue
        end do
    end do
    print*, dl(jml(2,2,1)), dl(jml(2,2,3)), it2
end subroutine find_boundary_layer_thickness


real(dp) function hyperviscosity(j)
    implicit none
    integer  :: j, jc
    real(dp) :: eta0, eta_max
    eta0    = 10.0_dp ** Ocean_viscosity
    eta_max = eta0 * 1e4_dp
    jc = max_jmax * 70/100
    if(SWITCH_hyperviscosity == 0)then
        hyperviscosity = 1.0_dp
    else
        if(j <= jc) then
            hyperviscosity = 1.0_dp
        else
            hyperviscosity = (eta_max/eta0) ** (real(j-jc, dp)/real(max_jmax - jc, dp))
        end if
    end if
end function hyperviscosity
end module mod_math