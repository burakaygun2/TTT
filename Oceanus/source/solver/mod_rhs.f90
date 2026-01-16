module mod_rhs
    use mod_parameters
    use mod_variables
    use mod_potential
    use mod_math
    use mod_SHTns
    implicit none
    
contains
    
subroutine generate_rhs(t2, jmind, sph, tor)
    implicit none
    integer    :: j, m, jm1, jp1, jm0, jms, smm, tmm, jmind
    integer    :: il, ir, nn, ilr
    real(dp)    :: dt_inertia, t2, jj, erjm1, erjp1, dtmu
    real(dp)    :: uc1, uc2, vc1, vc2, vc3, rg, dt2
    complex(dp) :: ur, avjm1, avjp1
    complex(dp) :: sph(:), tor(:)

    j = jmindx(jmind,2)
    m = jmindx(jmind,3)
    jm1 = jml(j,m,j-1)
    jm0 = jml(j,m,j  )
    jp1 = jml(j,m,j+1)
    jms =  jm(j,m)-1

    jj = real(j, dp)
    erjm1 =  sqrt(  jj           / (2.0_dp * jj + 1.0_dp))
    erjp1 = -sqrt( (jj+1.0_dp)   / (2.0_dp * jj + 1.0_dp))
    uc1 = -erjm1; uc2 = -erjp1;
    vc1 =       jj             / (2.0_dp * jj + 1.0_dp)
    vc2 =          (jj+1.0_dp)  / (2.0_dp * jj + 1.0_dp)
    vc3 = sqrt(jj*(jj+1.0_dp)) / (2.0_dp * jj + 1.0_dp)
    dt_inertia = -SWITCH_inertia/time_step
    dt2        = 0.5_dp * time_step
    smm = 0; tmm =0;
    sph(:) = dzero
    tor(:) = dzero

    if(abs(CN_theta - 1.0_dp)>1d-5) call Divergence_tensor
    ilr = 1
    do il = 1, numoflayers
        nn   = layer(il)%n
        do ir = 1, nn
            sph(smm+7*(ir-1)+1) = layer(il)%rhoe(ir, jms+1)/time_step
            ! layer(il)%rhoe(ir, jms+1)/time_step
            if(SWITCH_compressible==1) then
                sph(smm+7*(ir-1)+1) = (1 - layer_ocean_id(il))*layer(il)%rhoe(ir, jms+1)/time_step
                sph(smm+7*(ir-1)+2) = (1 - layer_ocean_id(il))*pjmstr*layer(il)%str (6*jms+1, ir)/(layer(il)%K_i(ir)*time_step)
            end if
            dtmu = 1.0_dp/(2.0_dp*layer(il)%mu_i(ir)*time_step)
            sph(smm+7*(ir-1)+5) = layer(il)%str (6*jms+2, ir)*dtmu ! !DEVIATORIC PART OF THE STRESS TENSOR TIME DERIVATIVE
            sph(smm+7*(ir-1)+6) = layer(il)%str (6*jms+4, ir)*dtmu ! !DEVIATORIC PART OF THE STRESS TENSOR TIME DERIVATIVE
            sph(smm+7*(ir-1)+7) = layer(il)%str (6*jms+6, ir)*dtmu ! !DEVIATORIC PART OF THE STRESS TENSOR TIME DERIVATIVE
            
            tor(tmm+3*(ir-1)+2) = layer(il)%str (6*jms+3, ir)*dtmu
            tor(tmm+3*(ir-1)+3) = layer(il)%str (6*jms+5, ir)*dtmu

            if(ir>1)then
                sph(smm+7*(ir-1)+3) = dt_inertia * layer(il)%vel(jm1, ir) * layer(il)%rho_l(ir)
                sph(smm+7*(ir-1)+4) = dt_inertia * layer(il)%vel(jp1, ir) * layer(il)%rho_l(ir)
                tor(tmm+3*(ir-1)+1) = dt_inertia * layer(il)%vel(jm0, ir) * layer(il)%rho_l(ir)
                if(SWITCH_self_grav == 1 )then
                    sph(smm+7*(ir-1)+3) = sph(smm+7*(ir-1)+3) - layer(il)%rho_l(ir) * layer(il)%grad_pot(jm1, ir-1)
                    sph(smm+7*(ir-1)+4) = sph(smm+7*(ir-1)+4) - layer(il)%rho_l(ir) * layer(il)%grad_pot(jp1, ir-1)
                end if
                if((SWITCH_coriolis == 1).and.(layer_ocean_id(il) == 1)) then
                    if (t2<2.5*time_step) then
                        sph(smm+7*(ir-1)+3) = sph(smm+7*(ir-1)+3) - (layer(il)%fc(jm1, ir))
                        sph(smm+7*(ir-1)+4) = sph(smm+7*(ir-1)+4) - (layer(il)%fc(jp1, ir))
                        tor(tmm+3*(ir-1)+1) = tor(tmm+3*(ir-1)+1) - (layer(il)%fc(jm0, ir))
                    else
                        sph(smm+7*(ir-1)+3) = sph(smm+7*(ir-1)+3) - (1.5_dp * layer(il)%fc(jm1, ir) - 0.5_dp * layer(il)%fc0(jm1, ir)) !* exp(-rotation_period/t2)
                        sph(smm+7*(ir-1)+4) = sph(smm+7*(ir-1)+4) - (1.5_dp * layer(il)%fc(jp1, ir) - 0.5_dp * layer(il)%fc0(jp1, ir)) !* exp(-rotation_period/t2)
                        tor(tmm+3*(ir-1)+1) = tor(tmm+3*(ir-1)+1) - (1.5_dp * layer(il)%fc(jm0, ir) - 0.5_dp * layer(il)%fc0(jm0, ir)) !* exp(-rotation_period/t2)
                    end if
                    layer(il)%fc0(jm1, ir) =  layer(il)%fc(jm1, ir)
                    layer(il)%fc0(jp1, ir) =  layer(il)%fc(jp1, ir)
                    layer(il)%fc0(jm0, ir) =  layer(il)%fc(jm0, ir)
                end if

                if((SWITCH_vgradv == 1).and.(1 == layer_ocean(il))) then
                    if (t2<2.5*time_step) then
                        sph(smm+7*(ir-1)+3) = sph  (smm+7*(ir-1)+3) + layer(il)%vgradv(ir-1, jm1)
                        sph(smm+7*(ir-1)+4) = sph  (smm+7*(ir-1)+4) + layer(il)%vgradv(ir-1, jp1)
                        tor(tmm+3*(ir-1)+1) = tor  (tmm+3*(ir-1)+1) + layer(il)%vgradv(ir-1, jm0)
                    else
                        sph(smm+7*(ir-1)+3) = sph  (smm+7*(ir-1)+3) + (1.5_dp * layer(il)%vgradv(ir-1, jm1) - 0.5_dp * layer(il)%vgradv0(ir-1, jm1))
                        sph(smm+7*(ir-1)+4) = sph  (smm+7*(ir-1)+4) + (1.5_dp * layer(il)%vgradv(ir-1, jp1) - 0.5_dp * layer(il)%vgradv0(ir-1, jm1))
                        tor(tmm+3*(ir-1)+1) = tor  (tmm+3*(ir-1)+1) + (1.5_dp * layer(il)%vgradv(ir-1, jm0) - 0.5_dp * layer(il)%vgradv0(ir-1, jm1))
                    end if
                    layer(il)%vgradv0(ir-1, jm1) =  layer(il)%vgradv(ir-1, jm1)
                    layer(il)%vgradv0(ir-1, jp1) =  layer(il)%vgradv(ir-1, jp1)
                    layer(il)%vgradv0(ir-1, jm0) =  layer(il)%vgradv(ir-1, jm0)
                end if
                if(j==2)then
                    if(SWITCH_mmtides==0) then
                        sph(smm+7*(ir-1)+3) = sph(smm+7*(ir-1)+3) + CN_theta * tides(layer(il)%rho_l(ir), layer(il)%r_l(ir), m, t2, 'g') + (1.0_dp - CN_theta) * tides(layer(il)%rho_l(ir), layer(il)%r_l(ir), m,  t2-time_step, 'g')
                    end if
                    if(SWITCH_mmtides==1) sph(smm+7*(ir-1)+3) = sph(smm+7*(ir-1)+3) + CN_theta * mm_tides(layer(il)%rho_l(ir), layer(il)%r_l(ir), m, t2, 'g') + (1.0_dp - CN_theta) * mm_tides(layer(il)%rho_l(ir), layer(il)%r_l(ir), m, t2-time_step, 'g')
                    if(SWITCH_mmtides==2) sph(smm+7*(ir-1)+3) = sph(smm+7*(ir-1)+3) + CN_theta * lunar_tides(layer(il)%rho_l(ir), layer(il)%r_l(ir), m, t2, 'g') + (1.0_dp - CN_theta) * lunar_tides(layer(il)%rho_l(ir), layer(il)%r_l(ir), m, t2-time_step, 'g')
                end if
            end if
        end do
        smm  = sum(layer(1:il)%sm); tmm = sum(layer(1:il)%tm)
    end do

    do il = 1, numoflayers-1
        smm   = sum(layer(1:il)%sm); tmm = sum(layer(1:il)%tm)
        nn    = layer(il)%n
        ur  = layer(il)%ur_i(jm(j,m), 2)
        
        avjm1 =     (layer(il)%vel(jm1, nn  ) * (layer(il)%r_l(nn+1) - layer(il)%r_i(nn)) &
                +    layer(il)%vel(jm1, nn+1) * (layer(il)%r_i(nn)   - layer(il)%r_l(nn)))/ &
                    (layer(il)%r_l(nn+1) - layer(il)%r_l(nn))
        avjp1 =  (layer(il)%vel(jp1, nn  ) * (layer(il)%r_l(nn+1) - layer(il)%r_i(nn)) &
                + layer(il)%vel(jp1, nn+1) * (layer(il)%r_i(nn)   - layer(il)%r_l(nn)))/ &
                (layer(il)%r_l(nn+1) - layer(il)%r_l(nn))

        ! !INTERFACE CONDITIONS (INTEGRATION OF THE RADIAL VELOCITY)
        rg = (layer(il)%rho_i(nn) - layer(il+1)%rho_i(1)) * layer(il)%g_i(nn)
        sph(smm+3) = (-rg * (uc1 * ur - dt2 * vc1 * avjm1 + dt2 * vc3 * avjp1))
        sph(smm+4) = (-rg * (uc2 * ur + dt2 * vc3 * avjm1 - dt2 * vc2 * avjp1))
    end do
    il = numoflayers; nn = layer(il)%n
    ur = layer(il)%ur_i(jm(j,m), 2)
    avjm1 = 0.5_dp * (layer(il)%vel(jm1, nn) + layer(il)%vel(jm1, nn+1))
    avjp1 = 0.5_dp * (layer(il)%vel(jp1, nn) + layer(il)%vel(jp1, nn+1))
    rg = (layer(numoflayers)%rho_i(nn)) * layer(numoflayers)%g_i(nn)
    if (SWITCH_one_layer == 'm') then
        sph(smtot-1) = ( rg * (uc1 * ur - dt2 * vc1 * avjm1 + dt2 * vc3 * avjp1))
        sph(smtot  ) = ( rg * (uc2 * ur + dt2 * vc3 * avjm1 - dt2 * vc2 * avjp1))
    else
        if (SWITCH_top_boundary_cond == 'fsv') then
            if((j==2).and.(m==2)) then
                ! print*, orb
                sph(smtot  ) = 3.86e-3_dp * sqrt(8.0_dp * pi / 15.0_dp) * (cos(orb * t2) - Im * sin(orb * t2))
            end if
        end if
    end if

end subroutine generate_rhs
end module mod_rhs