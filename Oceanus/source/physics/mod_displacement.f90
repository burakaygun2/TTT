module mod_displacement
    use mod_parameters
    use mod_variables
    implicit none
    
contains

subroutine boundary_displacement
    implicit none
    integer  :: j,  il, jm1, jp1, k, jmind, m
    real(dp) :: dt2, erjm1, erjp1, etjm1, etjp1, v_lower, v_upper
    complex(dp) :: avjm1, avjp1
    dt2 = time_step * 0.5_dp
    !!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(il,k,v_upper,v_lower,jmind,j,m,jm1,jp1,erjm1,erjp1,avjm1,avjp1) SCHEDULE(static)
    do il = 1, numoflayers
        k = layer(il)%n
        v_upper = (layer(il)%r_l(k+1) - layer(il)%r_i(k)) / (layer(il)%r_l(k+1) - layer(il)%r_l(k))
        v_lower = (layer(il)%r_i(k)   - layer(il)%r_l(k)) / (layer(il)%r_l(k+1) - layer(il)%r_l(k))
        do jmind = 2, compute_scalar_harm
            j = jmindx(jmind, 2); m   = jmindx(jmind, 3)
            jm1 = jml(j, m, j-1); jp1 = jml(j, m, j+1)
            erjm1 =  sqrt( real(j,dp)         / (2.0_dp*real(j,dp)+1.0_dp))
            erjp1 = -sqrt((real(j,dp)+1.0_dp) / (2.0_dp*real(j,dp)+1.0_dp))
            etjm1 =  sqrt((real(j,dp)+1.0_dp) / (2.0_dp*real(j,dp)+1.0_dp))
            etjp1 =  sqrt( real(j,dp)         / (2.0_dp*real(j,dp)+1.0_dp))
            avjm1 = v_upper*layer(il)%vel(jm1, k)+v_lower*layer(il)%vel(jm1, k+1)
            avjp1 = v_upper*layer(il)%vel(jp1, k)+v_lower*layer(il)%vel(jp1, k+1)

            layer(il)%ur_i(jmind, 2) = layer(il)%ur_i(jmind, 2) + (avjm1*erjm1 + avjp1*erjp1 + layer(il)%pvel (jmind, 2)) * dt2
            layer(il)%ut_i(jmind, 2) = layer(il)%ut_i(jmind, 2) + (avjm1*etjm1 + avjp1*etjp1 + layer(il)%pvelt(jmind, 2)) * dt2
            layer(il)%pvel (jmind, 2) = avjm1*erjm1 + avjp1*erjp1
            layer(il)%pvelt(jmind, 2) = avjm1*etjm1 + avjp1*etjp1
        end do
    end do
    !!$OMP END PARALLEL DO

end subroutine boundary_displacement


subroutine displacement
    implicit none
    integer :: j, m, i, jm1, jp1, jm0, k
    real(dp) :: v_upper, v_lower, dt2, erjm1, erjp1
    complex(dp) :: avjm1, avjp1, avj
    complex(dp) :: pavjm1, pavjp1, pavj

    dt2 = 0.5_dp * time_step
    do i = 1, numoflayers
        do j = jstart, jfinal
            erjm1 =  sqrt( real(j,dp)       / (2.0_dp*real(j,dp)+1.0_dp))
            erjp1 = -sqrt((real(j,dp)+1.0_dp)/ (2.0_dp*real(j,dp)+1.0_dp))
            do m = 0, min(j,2), 2
                jm1 = jml2(j,m,j-1)
                jm0 = jml2(j,m,j  )
                jp1 = jml2(j,m,j+1)
                do k = 1, layer(i)%n
                    v_upper = (layer(i)%r_l(k+1) - layer(i)%r_i(k)) / (layer(i)%r_l(k+1) - layer(i)%r_l(k))
                    v_lower = (layer(i)%r_i(k)   - layer(i)%r_l(k)) / (layer(i)%r_l(k+1) - layer(i)%r_l(k))
                    avjm1 = v_upper*layer(i)%vel(k, jm1)+v_lower*layer(i)%vel(k+1, jm1)
                    avj   = v_upper*layer(i)%vel(k, jm0)+v_lower*layer(i)%vel(k+1, jm0)
                    avjp1 = v_upper*layer(i)%vel(k, jp1)+v_lower*layer(i)%vel(k+1, jp1)

                    pavjm1 = v_upper*layer(i)%pvel(k, jm1)+v_lower*layer(i)%pvel(k+1, jm1)
                    pavj   = v_upper*layer(i)%pvel(k, jm0)+v_lower*layer(i)%pvel(k+1, jm0)
                    pavjp1 = v_upper*layer(i)%pvel(k, jp1)+v_lower*layer(i)%pvel(k+1, jp1)

                    layer(i)%dis(k, jm1) = layer(i)%dis(k, jm1) + dt2 * (avjm1 + pavjm1)
                    layer(i)%dis(k, jm0) = layer(i)%dis(k, jm0) + dt2 * (avj   + pavj  )
                    layer(i)%dis(k, jp1) = layer(i)%dis(k, jp1) + dt2 * (avjp1 + pavjp1)

                    layer(i)%ur_i(k, jm2(j,m)) = layer(i)%dis(k, jm1)*erjm1 + layer(i)%dis(k, jp1)*erjp1
                end do
            end do
        end do
    end do
end subroutine displacement

subroutine displacement2
    implicit none
    integer    :: j, m, il, jm1, jp1, jm0, ir
    real(dp)    :: v_upper, v_lower, dt2, erjm1, erjp1
    complex(dp) :: avjm1, avjp1, avj

    dt2 = 0.5_dp * time_step
    do il = 1, numoflayers
        do j = jstart, jfinal
            erjm1 =  sqrt( real(j,dp)       / (2.0_dp*real(j,dp)+1.0_dp))
            erjp1 = -sqrt((real(j,dp)+1.0_dp)/ (2.0_dp*real(j,dp)+1.0_dp))
            do m = 0, min(j,2), 2
                jm1 = jml2(j,m,j-1)
                jm0 = jml2(j,m,j  )
                jp1 = jml2(j,m,j+1)
                do ir = 1, layer(il)%n
                    v_upper = (layer(il)%r_l(ir+1) - layer(il)%r_i(ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
                    v_lower = (layer(il)%r_i(ir)   - layer(il)%r_l(ir)) / (layer(il)%r_l(ir+1) - layer(il)%r_l(ir))
                    avjm1 = v_upper*layer(il)%vel(ir, jm1)+v_lower*layer(il)%vel(ir+1, jm1)
                    avj   = v_upper*layer(il)%vel(ir, jm0)+v_lower*layer(il)%vel(ir+1, jm0)
                    avjp1 = v_upper*layer(il)%vel(ir, jp1)+v_lower*layer(il)%vel(ir+1, jp1)

                    layer(il)%dis(ir, jm1) = layer(il)%dis(ir, jm1) + dt2 * (avjm1)
                    layer(il)%dis(ir, jm0) = layer(il)%dis(ir, jm0) + dt2 * (avj  )
                    layer(il)%dis(ir, jp1) = layer(il)%dis(ir, jp1) + dt2 * (avjp1)

                    layer(il)%ur_i(ir, jm2(j,m)) = layer(il)%dis(ir, jm1)*erjm1 + layer(il)%dis(ir, jp1)*erjp1
                end do
            end do
        end do
    end do
end subroutine displacement2

!! for the first step of rerun only
subroutine radial_displacement
    implicit none
    integer :: j, m, i, jm1, jp1, k
    real(dp) :: erjm1, erjp1
    

    do i = 1, numoflayers
        do j = jstart, jfinal
            erjm1 =  sqrt( real(j,dp)       / (2.0_dp*real(j,dp)+1.0_dp))
            erjp1 = -sqrt((real(j,dp)+1.0_dp)/ (2.0_dp*real(j,dp)+1.0_dp))
            do m = 0, min(j,2), 2
                jm1 = jml2(j,m,j-1)
                jp1 = jml2(j,m,j+1)
                do k = 1, layer(i)%n
                    layer(i)%ur_i(k, jm2(j,m)) = layer(i)%dis(k, jm1)*erjm1 + layer(i)%dis(k, jp1)*erjp1
                end do
                
            end do
        end do
    end do
end subroutine radial_displacement

subroutine layer_displacement
    implicit none
    integer :: j, m, i, jm1, jp1, jm0, k
    real(dp) :: dt2, erjm1, erjp1
    complex(dp) :: avjm1, avjp1, avj
    complex(dp) :: pavjm1, pavjp1, pavj

    dt2 = 0.5_dp * time_step
    do i = 1, numoflayers
        do j = jstart, jfinal
            erjm1 =  sqrt( real(j,dp)       / (2.0_dp*real(j,dp)+1.0_dp))
            erjp1 = -sqrt((real(j,dp)+1.0_dp)/ (2.0_dp*real(j,dp)+1.0_dp))
            do m = 0, min(j,2), 2
                jm1 = jml2(j,m,j-1)
                jm0 = jml2(j,m,j  )
                jp1 = jml2(j,m,j+1)
                do k = 1, layer(i)%n+1
                    avjm1 = layer(i)%vel(k, jm1)
                    avj   = layer(i)%vel(k, jm0)
                    avjp1 = layer(i)%vel(k, jp1)

                    pavjm1 = layer(i)%pvel(k, jm1)
                    pavj   = layer(i)%pvel(k, jm0)
                    pavjp1 = layer(i)%pvel(k, jp1)

                    layer(i)%disl(k, jm1) = layer(i)%disl(k, jm1) + dt2 * (avjm1 + pavjm1)
                    layer(i)%disl(k, jm0) = layer(i)%disl(k, jm0) + dt2 * (avj   + pavj  )
                    layer(i)%disl(k, jp1) = layer(i)%disl(k, jp1) + dt2 * (avjp1 + pavjp1)
                end do
                
            end do
        end do
    end do
end subroutine layer_displacement

!!displacement formulation only:
subroutine rad_u_displacement
    implicit none
    integer :: j, m, i, jm1, jp1, jm0, k
    real(dp) :: v_upper, v_lower,  erjm1, erjp1
    complex(dp) :: avjm1, avjp1, avj 

    do i = 1, numoflayers
        do j = jstart, jfinal
            erjm1 =  sqrt( real(j,dp)       / (2.0_dp*real(j,dp)+1.0_dp))
            erjp1 = -sqrt((real(j,dp)+1.0_dp)/ (2.0_dp*real(j,dp)+1.0_dp))
            do m = 0, min(j,2), 2
                jm1 = jml2(j,m,j-1)
                jm0 = jml2(j,m,j  )
                jp1 = jml2(j,m,j+1)
                do k = 1, layer(i)%n
                    v_upper = (layer(i)%r_l(k+1) - layer(i)%r_i(k))  / (layer(i)%r_l(k+1) - layer(i)%r_l(k))
                    v_lower = (layer(i)%r_i(k)    - layer(i)%r_l(k)) / (layer(i)%r_l(k+1) - layer(i)%r_l(k))
                    avjm1 = v_upper*layer(i)%vel(k, jm1)+v_lower*layer(i)%vel(k+1, jm1)
                    avj   = v_upper*layer(i)%vel(k, jm0)+v_lower*layer(i)%vel(k+1, jm0)
                    avjp1 = v_upper*layer(i)%vel(k, jp1)+v_lower*layer(i)%vel(k+1, jp1)

                    layer(i)%ur_i(k, jm2(j,m)) = avjm1*erjm1 + avjp1*erjp1
                end do
                
            end do
        end do
    end do
end subroutine rad_u_displacement

subroutine update_maxwell_term
    implicit none
    integer :: il, ir, j, m, jjm2, jj2, jjp2


    do il = 1, numoflayers
        do j = jstart, jfinal
            do m = 0, min(j,2), 2
                jjm2 = 6*(jm2(j,m)-1)+2
                jj2  = 6*(jm2(j,m)-1)+4
                jjp2 = 6*(jm2(j,m)-1)+6
                do ir = 1, layer(il)%n
                    layer(il)%maxwell_term(ir, jjm2) = layer(il)%maxwell_term(ir, jjm2) + layer(il)%str(ir, jjm2)
                    layer(il)%maxwell_term(ir, jj2 ) = layer(il)%maxwell_term(ir, jj2 ) + layer(il)%str(ir, jj2 )
                    layer(il)%maxwell_term(ir, jjp2) = layer(il)%maxwell_term(ir, jjp2) + layer(il)%str(ir, jjp2)
                end do
            end do
        end do
    end do
    
end subroutine update_maxwell_term

subroutine update_andrade_term(tt, itime)
    implicit none
    real(dp) :: tt, weight1, weight2, dt2
    integer :: itime, il, ir, itt, j, m, jjm2, jj2, jjp2
    ! open(500, file='andrade_stress.dat', form='unformatted',access='stream',status='replace')


    dt2 = 0.5_dp * time_step
    if(itime == 1)then
        do il = 1, numoflayers
            do ir = 1, layer(il)%n
                layer(il)%andrade_CC(ir) = ((zeta*layer(il)%eta_i(ir))**(-alpha)) * (layer(il)%mu_i(ir)**(alpha-1))
            end do
        end do
    end if
    do il = 1, numoflayers
        layer(il)%andrade_term = dzero
    end do
    do itt = 1, itime-1
        if(itt == 1)then
            weight1 = 0.0_dp
        else
            weight1 = alpha * (tt -  itt   *time_step)**(alpha - 1)
        end if
        if((itt+1)==itime)then
            weight2 = 0.0_dp
        else
            weight2 = alpha * (tt - (itt+1)*time_step)**(alpha - 1)
        end if
        do il = 1, numoflayers
            do j = jstart, jfinal
                do m = 0, min(2,j),2
                    jjm2 = 6*(jm2(j,m)-1)+2
                    jj2  = 6*(jm2(j,m)-1)+4
                    jjp2 = 6*(jm2(j,m)-1)+6
                    layer(il)%andrade_term(:,jjm2) = layer(il)%andrade_term(:,jjm2) + dt2 * layer(il)%andrade_CC(:) * (weight1*layer(il)%str_a(:,jjm2,itt) + weight2*layer(il)%str_a(:,jjm2,itt+1))
                    layer(il)%andrade_term(:,jj2 ) = layer(il)%andrade_term(:,jj2 ) + dt2 * layer(il)%andrade_CC(:) * (weight1*layer(il)%str_a(:,jj2 ,itt) + weight2*layer(il)%str_a(:,jj2 ,itt+1))
                    layer(il)%andrade_term(:,jjp2) = layer(il)%andrade_term(:,jjp2) + dt2 * layer(il)%andrade_CC(:) * (weight1*layer(il)%str_a(:,jjp2,itt) + weight2*layer(il)%str_a(:,jjp2,itt+1))
                end do
            end do
        end do
    end do
end subroutine update_andrade_term
    
end module mod_displacement