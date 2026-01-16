module mod_matrix
    use mod_allocation
    use mod_variables
    use mod_parameters
    use mod_matrix_coeffs
    use mod_LU
    use mod_math
    implicit none
    
contains

subroutine create_matrix(j)
    implicit none
    integer :: j
    real(dp) :: jj, d,d2


    jj = real(j, dp)
    call matrix_coeffs(jj)
    S_U(:,:,j)=zero; T_U(:,:,j)=zero; S_L(:,:,j)=zero; T_L(:,:,j)=zero; S_I(:,j)=int(zero); T_I(:,j)=int(zero)
    if(SWITCH_compressible==0)then
        call generate_matrix_incompressible(j)
        call internal_boundaries(j)
        call top_bottom_boundaries(j)
    elseif(SWITCH_compressible==1)then
        call generate_matrix_compressible(j)
        call internal_boundaries(j)
        call top_bottom_boundaries(j)
    else
        stop "For incompressible: 0 or Compressible: 1..."
    end if
    call lu_decomp(S_U(:,:,j), smtot, sph_lo_diag, sph_hi_diag, smtot, sph_band, S_L(:,:,j), sph_lo_diag, S_I(:, j), d)
    call lu_decomp(T_U(:,:,j), tmtot, tor_lo_diag, tor_hi_diag, tmtot, tor_band, T_L(:,:,j), tor_lo_diag, T_I(:, j), d2)

end subroutine create_matrix
    
subroutine generate_matrix_incompressible(j)
    implicit none
    integer :: ir, il, smm, nn, tmm, j
    real(dp) :: hl, hup, hlo, hup_r, hlo_r
    real(dp) :: hi, hup2, hlo2, hup2_r, hlo2_r

    smm = 0; tmm=0
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn
            hl  =  layer(il)%r_l(ir+1) - layer(il)%r_l(ir)
            hup = (layer(il)%r_l(ir+1) - layer(il)%r_i(ir)) / hl;           hup_r = hup/layer(il)%r_i(ir)
            hlo = (layer(il)%r_i(ir)   - layer(il)%r_l(ir)) / hl;           hlo_r = hlo/layer(il)%r_i(ir)
            !! Balance of mass
            ! if (constant_density.eqv.(.TRUE.)) then
            !     S_U(smm+7*(ir-1)+1, sph_diag+2, j) = 1.0_dp                                        ! !RHO_t
            ! else

            S_U(smm+7*(ir-1)+1, sph_diag+2, j) = 1.0_dp/time_step                              ! !RHO_t
            S_U(smm+7*(ir-1)+1, sph_diag+0, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjm1*hup   ! !V_{JM}^{J-1}, LOWER
            S_U(smm+7*(ir-1)+1, sph_diag+1, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjp1*hup   ! !V_{JM}^{J+1}, LOWER
            S_U(smm+7*(ir-1)+1, sph_diag+7, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjm1*hlo   ! !V_{JM}^{J-1}, UPPER
            S_U(smm+7*(ir-1)+1, sph_diag+8, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjp1*hlo   ! !V_{JM}^{J+1}, UPPER
            ! end if
            !! Constitutive - pressure div(v) = 0
            S_U(smm+7*(ir-1)+2, sph_diag-1, j) = -sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hup_r
            S_U(smm+7*(ir-1)+2, sph_diag+0, j) = -sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hup_r
            S_U(smm+7*(ir-1)+2, sph_diag+6, j) =  sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hlo_r
            S_U(smm+7*(ir-1)+2, sph_diag+7, j) =  sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hlo_r
            !! CONSTITUTIVE EQUATION Y_{JM}^{J-2,2} 
            S_U(smm+7*(ir-1)+5, sph_diag-4, j) = (-sc%Const_jm2_vjm1(1)/hl + sc%Const_jm2_vjm1(2)*hup_r)
            S_U(smm+7*(ir-1)+5, sph_diag+0, j) = 1.0_dp/(2.0_dp*layer(il)%eta_i(ir)) + 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) 
            S_U(smm+7*(ir-1)+5, sph_diag+3, j) = ( sc%Const_jm2_vjm1(1)/hl + sc%Const_jm2_vjm1(2)*hlo_r)
            !! CONSTITUTIVE EQUATION Y_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+6, sph_diag-5, j) = (-sc%Const_j2__vjm1(1)/hl + sc%Const_j2__vjm1(2)*hup_r)
            S_U(smm+7*(ir-1)+6, sph_diag-4, j) = (-sc%Const_j2__vjp1(1)/hl + sc%Const_j2__vjp1(2)*hup_r)
            S_U(smm+7*(ir-1)+6, sph_diag+0, j) = 1.0_dp/(2.0_dp*layer(il)%eta_i(ir)) + 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) 
            S_U(smm+7*(ir-1)+6, sph_diag+2, j) = ( sc%Const_j2__vjm1(1)/hl + sc%Const_j2__vjm1(2)*hlo_r)
            S_U(smm+7*(ir-1)+6, sph_diag+3, j) = ( sc%Const_j2__vjp1(1)/hl + sc%Const_j2__vjp1(2)*hlo_r)
            !! CONSTITUTIVE EQUATION Y_{JM}^{J+2,2}
            S_U(smm+7*(ir-1)+7, sph_diag-5, j) = (-sc%Const_jp2_vjp1(1)/hl + sc%Const_jp2_vjp1(2)*hup_r)
            S_U(smm+7*(ir-1)+7, sph_diag+0, j) = 1.0_dp/(2.0_dp*layer(il)%eta_i(ir)) + 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) 
            S_U(smm+7*(ir-1)+7, sph_diag+2, j) = ( sc%Const_jp2_vjp1(1)/hl + sc%Const_jp2_vjp1(2)*hlo_r)
            if(SWITCH_ocean_status == 'v')then
                ! !CONSTITUTIVE EQUATION Y_{JM}^{J-1,2} 
                T_U(tmm+3*(ir-1) + 2, tor_diag-1, j) = (-tc%Const_jm1_vj(1)/hl + tc%Const_jm1_vj(2)*hup_r)
                T_U(tmm+3*(ir-1) + 2, tor_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + 1.0_dp/(2.0_dp*layer(il)%eta_i(ir))
                T_U(tmm+3*(ir-1) + 2, tor_diag+2, j) = ( tc%Const_jm1_vj(1)/hl + tc%Const_jm1_vj(2)*hlo_r)

                ! !CONSTITUTIVE EQUATION Y_{JM}^{J+1,2}
                T_U(tmm+3*(ir-1) + 3, tor_diag-2, j) = (-tc%Const_jp1_vj(1)/hl + tc%Const_jp1_vj(2)*hup_r)
                T_U(tmm+3*(ir-1) + 3, tor_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + 1.0_dp/(2.0_dp*layer(il)%eta_i(ir))
                T_U(tmm+3*(ir-1) + 3, tor_diag+1, j) = ( tc%Const_jp1_vj(1)/hl + tc%Const_jp1_vj(2)*hlo_r)
            end if
        end do
        smm = sum(layer(1:il)%sm); tmm=sum(layer(1:il)%tm)
    end do

    smm = 0; tmm=0
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 2, nn
            hi   =   layer(il)%r_i(ir)   - layer(il)%r_i(ir-1)
            hup2 =  (layer(il)%r_i(ir)   - layer(il)%r_l(  ir)) / hi;           hup2_r = hup2/layer(il)%r_l(ir)
            hlo2 =  (layer(il)%r_l(ir)   - layer(il)%r_i(ir-1)) / hi;           hlo2_r = hlo2/layer(il)%r_l(ir)
            !!Y_{JM}^{J-1} COMPONENT
            S_U(smm+7*(ir-1)+3, sph_diag-7, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jm1_rho          * hup2      ! !RHO_{JM}
            S_U(smm+7*(ir-1)+3, sph_diag-6, j) = (CN_theta)*(-sc%Mom__jm1_sj0(1) /hi + sc%Mom__jm1_sj0(2)  * hup2_r)   ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+3, sph_diag-5, j) = (CN_theta)*(-sc%Mom__jm1_sjm2(1)/hi + sc%Mom__jm1_sjm2(2) * hup2_r)   ! !S_{JM}^{J-2,2}
            S_U(smm+7*(ir-1)+3, sph_diag-4, j) = (CN_theta)*(-sc%Mom__jm1_sj2(1) /hi + sc%Mom__jm1_sj2(2)  * hup2_r)   ! !S_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+3, sph_diag+0, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jm1_rho          * hlo2      ! !RHO_{JM}
            S_U(smm+7*(ir-1)+3, sph_diag+1, j) = (CN_theta)*( sc%Mom__jm1_sj0(1) /hi + sc%Mom__jm1_sj0(2)  * hlo2_r)   ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+3, sph_diag+2, j) = (CN_theta)*( sc%Mom__jm1_sjm2(1)/hi + sc%Mom__jm1_sjm2(2) * hlo2_r)   ! !S_{JM}^{J-2,2}
            S_U(smm+7*(ir-1)+3, sph_diag+3, j) = (CN_theta)*( sc%Mom__jm1_sj2(1) /hi + sc%Mom__jm1_sj2(2)  * hlo2_r)   ! !S_{JM}^{J  ,2}

            S_U(smm+7*(ir-1)+3, sph_diag-2, j) =        -SWITCH_inertia*layer(il)%rho_l(ir)/time_step
            !!Y_{JM}^{J+1} COMPONENT
            S_U(smm+7*(ir-1)+4, sph_diag-8, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jp1_rho*hup2                ! !RHO_{JM}
            S_U(smm+7*(ir-1)+4, sph_diag-7, j) = (CN_theta)*(-sc%Mom__jp1_sj0(1) /hi + sc%Mom__jp1_sj0(2)  * hup2_r)  ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+4, sph_diag-5, j) = (CN_theta)*(-sc%Mom__jp1_sj2(1) /hi + sc%Mom__jp1_sj2(2)  * hup2_r)  ! !S_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+4, sph_diag-4, j) = (CN_theta)*(-sc%Mom__jp1_sjp2(1)/hi + sc%Mom__jp1_sjp2(2) * hup2_r)  ! !S_{JM}^{J+2,2}
            S_U(smm+7*(ir-1)+4, sph_diag-1, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jp1_rho*hlo2                ! !RHO_{JM}
            S_U(smm+7*(ir-1)+4, sph_diag+0, j) = (CN_theta)*( sc%Mom__jp1_sj0(1) /hi + sc%Mom__jp1_sj0(2)  * hlo2_r)  ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+4, sph_diag+2, j) = (CN_theta)*( sc%Mom__jp1_sj2(1) /hi + sc%Mom__jp1_sj2(2)  * hlo2_r)  ! !S_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+4, sph_diag+3, j) = (CN_theta)*( sc%Mom__jp1_sjp2(1)/hi + sc%Mom__jp1_sjp2(2) * hlo2_r)  ! !S_{JM}^{J+2,2}

            S_U(smm+7*(ir-1)+4, sph_diag-2, j) =        -SWITCH_inertia*layer(il)%rho_l(ir)/time_step

            if(SWITCH_ocean_status == 'v')then
                T_U(tmm+3*(ir-1)+1, tor_diag-2, j) = (CN_theta)*(-tc%Mom__j_s_jm1(1) /hi + tc%Mom__j_s_jm1(2)  * hup2_r)
                T_U(tmm+3*(ir-1)+1, tor_diag-1, j) = (CN_theta)*(-tc%Mom__j_s_jp1(1) /hi + tc%Mom__j_s_jp1(2)  * hup2_r)
                T_U(tmm+3*(ir-1)+1, tor_diag+0, j) = -SWITCH_inertia*layer(il)%rho_l(ir)/time_step
                T_U(tmm+3*(ir-1)+1, tor_diag+1, j) = (CN_theta)*( tc%Mom__j_s_jm1(1) /hi + tc%Mom__j_s_jm1(2)  * hlo2_r)
                T_U(tmm+3*(ir-1)+1, tor_diag+2, j) = (CN_theta)*( tc%Mom__j_s_jp1(1) /hi + tc%Mom__j_s_jp1(2)  * hlo2_r)
            end if
        end do
        smm = sum(layer(1:il)%sm); tmm=sum(layer(1:il)%tm)
    end do
end subroutine generate_matrix_incompressible

subroutine generate_matrix_compressible(j)
    implicit none
    integer :: ir, il, smm, tmm, nn, j
    real(dp) :: hl, hup, hlo, hup_r, hlo_r
    real(dp) :: hi, hup2, hlo2, hup2_r, hlo2_r
    real(dp) :: anelas, vis
    
    smm = 0; tmm=0
    do il = 1, numoflayers
        nn = layer(il)%n
        if(j == 2) print*, layer_ocean_id(il), il
        if(layer_ocean_id(il) == 1)then
            anelas=0.0_dp
        else
            anelas=1.0_dp
        end if
        do ir = 1, nn
            hl  =  layer(il)%r_l(ir+1) - layer(il)%r_l(ir)
            hup = (layer(il)%r_l(ir+1) - layer(il)%r_i(ir)) / hl;           hup_r = hup/layer(il)%r_i(ir)
            hlo = (layer(il)%r_i(ir)   - layer(il)%r_l(ir)) / hl;           hlo_r = hlo/layer(il)%r_i(ir)
            !!Balance of mass
            S_U(smm+7*(ir-1)+1, sph_diag+0, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjm1*hup + layer(il)%rho_i(ir)*(-sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hup_r)  ! !V_{JM}^{J-1}, LOWER
            S_U(smm+7*(ir-1)+1, sph_diag+1, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjp1*hup + layer(il)%rho_i(ir)*(-sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hup_r)  ! !V_{JM}^{J+1}, LOWER
            S_U(smm+7*(ir-1)+1, sph_diag+2, j) = anelas/time_step                              ! !RHO_t
            S_U(smm+7*(ir-1)+1, sph_diag+7, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjm1*hlo + layer(il)%rho_i(ir)*( sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hlo_r)  ! !V_{JM}^{J-1}, UPPER
            S_U(smm+7*(ir-1)+1, sph_diag+8, j) = layer(il)%rho_r_i(ir)*sc%mass_bal_vjp1*hlo + layer(il)%rho_i(ir)*( sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hlo_r)  ! !V_{JM}^{J+1}, UPPER


            if(layer_ocean_id(il) == 1)then
                !! Constitutive - pressure p - K rho/rho0 + zeta div(v) = 0
                S_U(smm+7*(ir-1)+2, sph_diag+2, j) = pjmstr/(layer(il)%K_i(ir))   !! pressure
                S_U(smm+7*(ir-1)+2, sph_diag+1, j) = -1.0_dp /(layer(il)%rho_i(ir)) !! density
                S_U(smm+7*(ir-1)+2, sph_diag-1, j) = -0.0_dp*(-sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hup_r)
                S_U(smm+7*(ir-1)+2, sph_diag+0, j) = -0.0_dp*(-sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hup_r)
                S_U(smm+7*(ir-1)+2, sph_diag+6, j) = -0.0_dp*( sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hlo_r)
                S_U(smm+7*(ir-1)+2, sph_diag+7, j) = -0.0_dp*( sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hlo_r)
            else
                !! Constitutive - pressure p_t - vr rho0 g0 + K div(v) = 0
                S_U(smm+7*(ir-1)+2, sph_diag+2, j) = pjmstr/(layer(il)%K_i(ir)*time_step)
                S_U(smm+7*(ir-1)+2, sph_diag-1, j) = (-sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hup_r) - (layer(il)%rho_i(ir)*layer(il)%g_i(ir)*sc%mass_bal_vjm1*hup) * (1.0_dp/layer(il)%K_i(ir))
                S_U(smm+7*(ir-1)+2, sph_diag+0, j) = (-sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hup_r) - (layer(il)%rho_i(ir)*layer(il)%g_i(ir)*sc%mass_bal_vjp1*hup) * (1.0_dp/layer(il)%K_i(ir))
                S_U(smm+7*(ir-1)+2, sph_diag+6, j) = ( sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hlo_r) - (layer(il)%rho_i(ir)*layer(il)%g_i(ir)*sc%mass_bal_vjm1*hlo) * (1.0_dp/layer(il)%K_i(ir))
                S_U(smm+7*(ir-1)+2, sph_diag+7, j) = ( sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hlo_r) - (layer(il)%rho_i(ir)*layer(il)%g_i(ir)*sc%mass_bal_vjp1*hlo) * (1.0_dp/layer(il)%K_i(ir))
            end if
            if(layer_ocean_id(il) == 1) then
                vis = 1.0_dp/(2.0_dp*layer(il)%eta_i(ir)*hyperviscosity(j))
            else
                vis = 1.0_dp/(2.0_dp*layer(il)%eta_i(ir))
            end if
            ! !CONSTITUTIVE EQUATION Y_{JM}^{J-2,2} 
            S_U(smm+7*(ir-1)+5, sph_diag-4, j) = (-sc%Const_jm2_vjm1(1)/hl + sc%Const_jm2_vjm1(2)*hup_r)
            S_U(smm+7*(ir-1)+5, sph_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + vis
            S_U(smm+7*(ir-1)+5, sph_diag+3, j) = ( sc%Const_jm2_vjm1(1)/hl + sc%Const_jm2_vjm1(2)*hlo_r)
            ! !CONSTITUTIVE EQUATION Y_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+6, sph_diag-5, j) = (-sc%Const_j2__vjm1(1)/hl + sc%Const_j2__vjm1(2)*hup_r)
            S_U(smm+7*(ir-1)+6, sph_diag-4, j) = (-sc%Const_j2__vjp1(1)/hl + sc%Const_j2__vjp1(2)*hup_r)
            S_U(smm+7*(ir-1)+6, sph_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + vis
            S_U(smm+7*(ir-1)+6, sph_diag+2, j) = ( sc%Const_j2__vjm1(1)/hl + sc%Const_j2__vjm1(2)*hlo_r)
            S_U(smm+7*(ir-1)+6, sph_diag+3, j) = ( sc%Const_j2__vjp1(1)/hl + sc%Const_j2__vjp1(2)*hlo_r)
            ! !CONSTITUTIVE EQUATION Y_{JM}^{J+2,2}
            S_U(smm+7*(ir-1)+7, sph_diag-5, j) = (-sc%Const_jp2_vjp1(1)/hl + sc%Const_jp2_vjp1(2)*hup_r)
            S_U(smm+7*(ir-1)+7, sph_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + vis
            S_U(smm+7*(ir-1)+7, sph_diag+2, j) = ( sc%Const_jp2_vjp1(1)/hl + sc%Const_jp2_vjp1(2)*hlo_r)
            if(SWITCH_ocean_status == 'v')then
                ! !CONSTITUTIVE EQUATION Y_{JM}^{J-1,2} 
                T_U(tmm+3*(ir-1) + 2, tor_diag-1, j) = (-tc%Const_jm1_vj(1)/hl + tc%Const_jm1_vj(2)*hup_r)
                T_U(tmm+3*(ir-1) + 2, tor_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + vis
                T_U(tmm+3*(ir-1) + 2, tor_diag+2, j) = ( tc%Const_jm1_vj(1)/hl + tc%Const_jm1_vj(2)*hlo_r)

                ! !CONSTITUTIVE EQUATION Y_{JM}^{J+1,2}
                T_U(tmm+3*(ir-1) + 3, tor_diag-2, j) = (-tc%Const_jp1_vj(1)/hl + tc%Const_jp1_vj(2)*hup_r)
                T_U(tmm+3*(ir-1) + 3, tor_diag+0, j) = 1.0_dp / (2.0_dp * layer(il)%mu_i(ir)*time_step) + vis
                T_U(tmm+3*(ir-1) + 3, tor_diag+1, j) = ( tc%Const_jp1_vj(1)/hl + tc%Const_jp1_vj(2)*hlo_r)
            end if
        end do
        smm=sum(layer(1:il)%sm); tmm=sum(layer(1:il)%tm)
    end do

    smm = 0; tmm=0
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 2, nn
            hi   =   layer(il)%r_i(ir)   - layer(il)%r_i(ir-1)
            hup2 =  (layer(il)%r_i(ir)   - layer(il)%r_l(  ir)) / hi;           hup2_r = hup2/layer(il)%r_l(ir)
            hlo2 =  (layer(il)%r_l(ir)   - layer(il)%r_i(ir-1)) / hi;           hlo2_r = hlo2/layer(il)%r_l(ir)
            !!Y_{JM}^{J-1} COMPONENT
            S_U(smm+7*(ir-1)+3, sph_diag-7, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jm1_rho          * hup2      ! !RHO_{JM}
            S_U(smm+7*(ir-1)+3, sph_diag-6, j) = (CN_theta)*(-sc%Mom__jm1_sj0(1) /hi + sc%Mom__jm1_sj0(2)  * hup2_r)   ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+3, sph_diag-5, j) = (CN_theta)*(-sc%Mom__jm1_sjm2(1)/hi + sc%Mom__jm1_sjm2(2) * hup2_r)   ! !S_{JM}^{J-2,2}
            S_U(smm+7*(ir-1)+3, sph_diag-4, j) = (CN_theta)*(-sc%Mom__jm1_sj2(1) /hi + sc%Mom__jm1_sj2(2)  * hup2_r)   ! !S_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+3, sph_diag+0, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jm1_rho          * hlo2      ! !RHO_{JM}
            S_U(smm+7*(ir-1)+3, sph_diag+1, j) = (CN_theta)*( sc%Mom__jm1_sj0(1) /hi + sc%Mom__jm1_sj0(2)  * hlo2_r)   ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+3, sph_diag+2, j) = (CN_theta)*( sc%Mom__jm1_sjm2(1)/hi + sc%Mom__jm1_sjm2(2) * hlo2_r)   ! !S_{JM}^{J-2,2}
            S_U(smm+7*(ir-1)+3, sph_diag+3, j) = (CN_theta)*( sc%Mom__jm1_sj2(1) /hi + sc%Mom__jm1_sj2(2)  * hlo2_r)   ! !S_{JM}^{J  ,2}

            S_U(smm+7*(ir-1)+3, sph_diag-2, j) =        -SWITCH_inertia*layer(il)%rho_l(ir)/time_step
            !!Y_{JM}^{J+1} COMPONENT
            S_U(smm+7*(ir-1)+4, sph_diag-8, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jp1_rho*hup2                ! !RHO_{JM}
            S_U(smm+7*(ir-1)+4, sph_diag-7, j) = (CN_theta)*(-sc%Mom__jp1_sj0(1) /hi + sc%Mom__jp1_sj0(2)  * hup2_r)  ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+4, sph_diag-5, j) = (CN_theta)*(-sc%Mom__jp1_sj2(1) /hi + sc%Mom__jp1_sj2(2)  * hup2_r)  ! !S_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+4, sph_diag-4, j) = (CN_theta)*(-sc%Mom__jp1_sjp2(1)/hi + sc%Mom__jp1_sjp2(2) * hup2_r)  ! !S_{JM}^{J+2,2}
            S_U(smm+7*(ir-1)+4, sph_diag-1, j) = (CN_theta)*(layer(il)%g_l(ir)) * sc%Mom__jp1_rho*hlo2                ! !RHO_{JM}
            S_U(smm+7*(ir-1)+4, sph_diag+0, j) = (CN_theta)*( sc%Mom__jp1_sj0(1) /hi + sc%Mom__jp1_sj0(2)  * hlo2_r)  ! !S_{JM}^{J  ,0}
            S_U(smm+7*(ir-1)+4, sph_diag+2, j) = (CN_theta)*( sc%Mom__jp1_sj2(1) /hi + sc%Mom__jp1_sj2(2)  * hlo2_r)  ! !S_{JM}^{J  ,2}
            S_U(smm+7*(ir-1)+4, sph_diag+3, j) = (CN_theta)*( sc%Mom__jp1_sjp2(1)/hi + sc%Mom__jp1_sjp2(2) * hlo2_r)  ! !S_{JM}^{J+2,2}

            S_U(smm+7*(ir-1)+4, sph_diag-2, j) =        -SWITCH_inertia*layer(il)%rho_l(ir)/time_step

            if(SWITCH_ocean_status == 'v')then
                T_U(tmm+3*(ir-1)+1, tor_diag-2, j) = (CN_theta)*(-tc%Mom__j_s_jm1(1) /hi + tc%Mom__j_s_jm1(2)  * hup2_r)
                T_U(tmm+3*(ir-1)+1, tor_diag-1, j) = (CN_theta)*(-tc%Mom__j_s_jp1(1) /hi + tc%Mom__j_s_jp1(2)  * hup2_r)
                T_U(tmm+3*(ir-1)+1, tor_diag+0, j) = -SWITCH_inertia*layer(il)%rho_l(ir)/time_step
                T_U(tmm+3*(ir-1)+1, tor_diag+1, j) = (CN_theta)*( tc%Mom__j_s_jm1(1) /hi + tc%Mom__j_s_jm1(2)  * hlo2_r)
                T_U(tmm+3*(ir-1)+1, tor_diag+2, j) = (CN_theta)*( tc%Mom__j_s_jp1(1) /hi + tc%Mom__j_s_jp1(2)  * hlo2_r)
            end if
        end do
        smm = sum(layer(1:il)%sm); tmm=sum(layer(1:il)%tm)
    end do
end subroutine generate_matrix_compressible

subroutine internal_boundaries(j)
    implicit none
    integer :: il, nn, smm, tmm, j
    real(dp) :: hl, hup, hlo, hup_r, hlo_r, drho, g

    do il = 2, numoflayers
        nn = layer(il-1)%n; smm = sum(layer(1:il-1)%sm); tmm = sum(layer(1:il-1)%tm)
        drho  = (layer(il-1)%rho_i(layer(il-1)%n) - layer(il)%rho_i(1))
        g     = layer(il-1)%  g_i(layer(il-1)%n)
        hl    =  layer(il-1)%r_l(nn+1) - layer(il-1)%r_l(nn)
        hup   = (layer(il-1)%r_l(nn+1) - layer(il-1)%r_i(nn))/(hl);            hup_r  = hup/layer(il-1)%r_i(nn)
        hlo   = (layer(il-1)%r_i(nn)   - layer(il-1)%r_l(nn))/(hl);            hlo_r  = hlo/layer(il-1)%r_i(nn)
        ! !CONTIUNITY OF TRACTION_{JM}^{J-1}
        S_U(smm+3, sph_diag-11, j) =  -hup * sc_b%Bound_jm1_vjm1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J-1}, LOWER
        S_U(smm+3, sph_diag-10, j) =  -hup * sc_b%Bound_jm1_vjp1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J-1}, LOWER
        S_U(smm+3, sph_diag-8 , j) =  -sc_b%Bound_jm1_s_0   ! !S_{JM}^{J  ,0}, LOWER
        S_U(smm+3, sph_diag-7 , j) =  -sc_b%Bound_jm1_s_jm2 ! !S_{JM}^{J-2,2}, LOWER
        S_U(smm+3, sph_diag-6 , j) =  -sc_b%Bound_jm1_s_j2  ! !S_{JM}^{J  ,2}, LOWER
        S_U(smm+3, sph_diag-4 , j) =  -hlo * sc_b%Bound_jm1_vjm1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J-1}, LOWER
        S_U(smm+3, sph_diag-3 , j) =  -hlo * sc_b%Bound_jm1_vjp1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J-1}, LOWER
        S_U(smm+3, sph_diag+1 , j) =   sc_b%Bound_jm1_s_0   ! !S_{JM}^{J  ,0}, UPPER
        S_U(smm+3, sph_diag+2 , j) =   sc_b%Bound_jm1_s_jm2 ! !S_{JM}^{J-2,2}, UPPER
        S_U(smm+3, sph_diag+3 , j) =   sc_b%Bound_jm1_s_j2  ! !S_{JM}^{J  ,2}, UPPER
        ! ?==========================================================
        ! ?==========================================================
        ! !CONTIUNITY OF TRACTION_{JM}^{J-1}
        S_U(smm+4, sph_diag-12, j) =  -hup * sc_b%Bound_jp1_vjm1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J-1}, LOWER
        S_U(smm+4, sph_diag-11, j) =  -hup * sc_b%Bound_jp1_vjp1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J+1}, LOWER
        S_U(smm+4, sph_diag-9 , j) =  -sc_b%Bound_jp1_s_0   ! !S_{JM}^{J  ,0}, LOWER
        S_U(smm+4, sph_diag-7 , j) =  -sc_b%Bound_jp1_s_j2  ! !S_{JM}^{J-2,2}, LOWER
        S_U(smm+4, sph_diag-6 , j) =  -sc_b%Bound_jp1_s_jp2 ! !S_{JM}^{J  ,2}, LOWER
        S_U(smm+4, sph_diag-5 , j) =  -hlo * sc_b%Bound_jp1_vjm1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J-1}, LOWER
        S_U(smm+4, sph_diag-4 , j) =  -hlo * sc_b%Bound_jp1_vjp1 * drho * g * 0.5_dp * time_step ! !V_{JM}^{J+1}, LOWER
        S_U(smm+4, sph_diag+0 , j) =   sc_b%Bound_jp1_s_0   ! !S_{JM}^{J  ,0}, UPPER
        S_U(smm+4, sph_diag+2 , j) =   sc_b%Bound_jp1_s_j2  ! !S_{JM}^{J  ,2}, UPPER
        S_U(smm+4, sph_diag+3 , j) =   sc_b%Bound_jp1_s_jp2 ! !S_{JM}^{J+2,2}, UPPER

        T_U(tmm+1, tor_diag-3, j) = -tc_b%Bound_j_s_jm1 !core
        T_U(tmm+1, tor_diag-2, j) = -tc_b%Bound_j_s_jp1 !core
        T_U(tmm+1, tor_diag+1, j) =  tc_b%Bound_j_s_jm1!ocean
        T_U(tmm+1, tor_diag+2, j) =  tc_b%Bound_j_s_jp1!ocean
    end do
    do il = 1, numoflayers-1
        smm = sum(layer(1:il)%sm); tmm = sum(layer(1:il)%tm)
        ! !CONTIUNITY OF V_{JM}^{J-1}
        S_U(smm-1, sph_diag-7, j) = -0.5_dp ! !BELOW LAYER
        S_U(smm-1, sph_diag+0, j) = -0.5_dp ! !BELOW LAYER
        S_U(smm-1, sph_diag+2, j) =  0.5_dp ! !UPPER LAYER
        S_U(smm-1, sph_diag+9, j) =  0.5_dp ! !UPPER LAYER

        ! !CONTIUNITY OF V_{JM}^{J+1}
        S_U(smm  , sph_diag-7, j) = -0.5_dp ! !BELOW LAYER
        S_U(smm  , sph_diag+0, j) = -0.5_dp ! !BELOW LAYER
        S_U(smm  , sph_diag+2, j) =  0.5_dp ! !UPPER LAYER
        S_U(smm  , sph_diag+9, j) =  0.5_dp ! !UPPER LAYER
        if(SWITCH_ocean_status=='v')then
            T_U(tmm, tor_diag-3, j) = -0.5_dp
            T_U(tmm, tor_diag+0, j) = -0.5_dp
            T_U(tmm, tor_diag+1, j) =  0.5_dp
            T_U(tmm, tor_diag+4, j) =  0.5_dp
        end if
    end do
end subroutine internal_boundaries

subroutine top_bottom_boundaries(j)
    implicit none
    integer :: il, nn, j
    real(dp) :: hl, hup, hlo, rho, g, jj, cv1, cv2

    jj = real(j, dp)
    if(SWITCH_one_layer == 'm') then
        T_U(1, tor_diag+0, j) = 1.0_dp

        S_U(3, sph_diag-2, j) = 1.0_dp
        S_U(4, sph_diag-2, j) = 1.0_dp
    else
        if     (SWITCH_bottom_boundary_cond == 'fs' ) then
            !! Tangnetial traction zero
            S_U(3, sph_diag+2, j) = sc_b%Bound_sigmatan_jm2
            S_U(3, sph_diag+3, j) = sc_b%Bound_sigmatan_j2
            S_U(3, sph_diag+4, j) = sc_b%Bound_sigmatan_jp2
            cv1 =  sqrt(jj)
            cv2 = -sqrt(jj+1.0_dp)
            !! Radial velocity zero
            S_U(4, sph_diag-3, j) = 0.5_dp * cv1
            S_U(4, sph_diag-2, j) = 0.5_dp * cv2
            S_U(4, sph_diag+4, j) = 0.5_dp * cv1
            S_U(4, sph_diag+5, j) = 0.5_dp * cv2

            T_U(1, tor_diag+1, j) = tc_b%Bound_j_s_jm1
            T_U(1, tor_diag+2, j) = tc_b%Bound_j_s_jp1
        elseif (SWITCH_bottom_boundary_cond == 'fsf') then
            !! Tangnetial traction zero
            S_U(3, sph_diag+2, j) = sc_b%Bound_sigmatan_jm2
            S_U(3, sph_diag+3, j) = sc_b%Bound_sigmatan_j2
            S_U(3, sph_diag+4, j) = sc_b%Bound_sigmatan_jp2

            !! Radial traction is prescribed
            S_U(4, sph_diag+0, j) = -1.0_dp / sqrt(3.0_dp)
            S_U(4, sph_diag+1, j) = sc_b%Bound_sigmarad_jm2
            S_U(4, sph_diag+2, j) = sc_b%Bound_sigmarad_j2
            S_U(4, sph_diag+3, j) = sc_b%Bound_sigmarad_jp2

            T_U(1, tor_diag+1, j) = tc_b%Bound_j_s_jm1
            T_U(1, tor_diag+2, j) = tc_b%Bound_j_s_jp1
        elseif (SWITCH_bottom_boundary_cond == 'ns' ) then
            !! All velocity zero
            S_U(3, sph_diag-2, j) = 0.5_dp
            S_U(4, sph_diag-2, j) = 0.5_dp

            S_U(3, sph_diag+5, j) = 0.5_dp
            S_U(4, sph_diag+5, j) = 0.5_dp

            T_U(1, tor_diag+0, j) = 0.5_dp
            T_U(1, tor_diag+3, j) = 0.5_dp
        elseif (SWITCH_bottom_boundary_cond == 'nsf') then
            cv1 =            jj + 1.0_dp
            cv2 = sqrt(jj * (jj + 1.0_dp))
            !! Tangential velocity zero
            S_U(3, sph_diag-2, j) = 0.5_dp * cv1
            S_U(3, sph_diag-1, j) = 0.5_dp * cv2

            S_U(3, sph_diag+5, j) = 0.5_dp * cv1
            S_U(3, sph_diag+6, j) = 0.5_dp * cv2
            !! Radial traction is presbribed
            S_U(4, sph_diag+0, j) = -1.0_dp / sqrt(3.0_dp)
            S_U(4, sph_diag+1, j) = sc_b%Bound_sigmarad_jm2
            S_U(4, sph_diag+2, j) = sc_b%Bound_sigmarad_j2
            S_U(4, sph_diag+3, j) = sc_b%Bound_sigmarad_jp2

            T_U(1, tor_diag+0, j) = 0.5_dp
            T_U(1, tor_diag+3, j) = 0.5_dp
        else
            print*, 'Invalid boundary condition at the bottom boundary...'
            stop
        end if
    end if


    il = numoflayers
    nn = layer(il)%n
    rho =  layer(il)%rho_i(nn)
    g   =  layer(il)%g_i(nn)
    hl  =  layer(il)%r_l(nn+1) - layer(il)%r_l(nn)
    hup = (layer(il)%r_l(nn+1) - layer(il)%r_i(nn))/(hl)
    hlo = (layer(il)%r_i(nn)   - layer(il)%r_l(nn))/(hl)
    if (SWITCH_one_layer == 'm') then
        ! !Free surface: Te_r = - u_r * rho * g = - (u_i + (v_{i+1} + v_{i})*dt/2) * rho * g
        ! !Traction vector j-1
        S_U(smtot-1, sph_diag-7, j) = sc_b%Bound_jm1_vjm1 * hup * rho * g * 0.5_dp * time_step ! 
        S_U(smtot-1, sph_diag-6, j) = sc_b%Bound_jm1_vjp1 * hup * rho * g * 0.5_dp * time_step ! 
        S_U(smtot-1, sph_diag-4, j) = sc_b%Bound_jm1_s_0    ! !S_{JM}^{J  ,0}
        S_U(smtot-1, sph_diag-3, j) = sc_b%Bound_jm1_s_jm2  ! !S_{JM}^{J-2,2}
        S_U(smtot-1, sph_diag-2, j) = sc_b%Bound_jm1_s_j2   ! !S_{JM}^{J  ,2}
        S_U(smtot-1, sph_diag+0, j) = sc_b%Bound_jm1_vjm1* hlo * rho * g * 0.5_dp * time_step !
        S_U(smtot-1, sph_diag+1, j) = sc_b%Bound_jm1_vjp1* hlo * rho * g * 0.5_dp * time_step !


        ! !Traction vector j+1
        S_U(smtot  , sph_diag-8, j) = sc_b%Bound_jp1_vjm1* hup * rho * g * 0.5_dp * time_step !
        S_U(smtot  , sph_diag-7, j) = sc_b%Bound_jp1_vjp1* hup * rho * g * 0.5_dp * time_step !
        S_U(smtot  , sph_diag-5, j) = sc_b%Bound_jp1_s_0   ! !S_{JM}^{J  ,0}
        S_U(smtot  , sph_diag-3, j) = sc_b%Bound_jp1_s_j2  ! !S_{JM}^{J  ,2}
        S_U(smtot  , sph_diag-2, j) = sc_b%Bound_jp1_s_jp2 ! !S_{JM}^{J+2,2}
        S_U(smtot  , sph_diag-1, j) = sc_b%Bound_jp1_vjm1 * hlo * rho * g * 0.5_dp * time_step !
        S_U(smtot  , sph_diag+0, j) = sc_b%Bound_jp1_vjp1 * hlo * rho * g * 0.5_dp * time_step !


        T_U(tmtot, tor_diag-2, j) = tc_b%Bound_j_s_jm1
        T_U(tmtot, tor_diag-1, j) = tc_b%Bound_j_s_jp1
    else
        if     (SWITCH_top_boundary_cond == 'fs' ) then
            !! Tangnetial traction zero
            S_U(smtot-1, sph_diag-4, j) = sc_b%Bound_sigmatan_jm2
            S_U(smtot-1, sph_diag-3, j) = sc_b%Bound_sigmatan_j2
            S_U(smtot-1, sph_diag-2, j) = sc_b%Bound_sigmatan_jp2
            cv1 =  sqrt(jj)
            cv2 = -sqrt(jj+1.0_dp)
            !! Radial velocity zero
            S_U(smtot  , sph_diag-8, j) = 0.5_dp * cv1
            S_U(smtot  , sph_diag-7, j) = 0.5_dp * cv2
            S_U(smtot  , sph_diag-1, j) = 0.5_dp * cv1
            S_U(smtot  , sph_diag+0, j) = 0.5_dp * cv2

            T_U(tmtot, tor_diag-2, j) = tc_b%Bound_j_s_jm1
            T_U(tmtot, tor_diag-1, j) = tc_b%Bound_j_s_jp1
        elseif (SWITCH_top_boundary_cond == 'ns' ) then
            S_U(smtot-1, sph_diag-7, j) = 0.5_dp
            S_U(smtot-1, sph_diag+0, j) = 0.5_dp
            S_U(smtot  , sph_diag-7, j) = 0.5_dp
            S_U(smtot  , sph_diag+0, j) = 0.5_dp
        elseif (SWITCH_top_boundary_cond == 'fsf') then
            !! Tangnetial traction zero
            S_U(smtot-1, sph_diag-3, j) = sc_b%Bound_sigmatan_jm2
            S_U(smtot-1, sph_diag-2, j) = sc_b%Bound_sigmatan_j2
            S_U(smtot-1, sph_diag-1, j) = sc_b%Bound_sigmatan_jp2

            !! Radial traction is prescribed
            S_U(smtot  , sph_diag-5, j) = -1.0_dp / sqrt(3.0_dp)
            S_U(smtot  , sph_diag-4, j) = sc_b%Bound_sigmarad_jm2
            S_U(smtot  , sph_diag-3, j) = sc_b%Bound_sigmarad_j2
            S_U(smtot  , sph_diag-2, j) = sc_b%Bound_sigmarad_jp2

            T_U(tmtot, tor_diag-2, j) = tc_b%Bound_j_s_jm1
            T_U(tmtot, tor_diag-1, j) = tc_b%Bound_j_s_jp1
        elseif (SWITCH_top_boundary_cond == 'nsf') then
            cv1 =            jj + 1.0_dp
            cv2 = sqrt(jj * (jj + 1.0_dp))
            !! Tangential velocity zero
            S_U(smtot-1, sph_diag+0, j) = 0.5_dp * cv1
            S_U(smtot-1, sph_diag+1, j) = 0.5_dp * cv2
            S_U(smtot-1, sph_diag-7, j) = 0.5_dp * cv1
            S_U(smtot-1, sph_diag-6, j) = 0.5_dp * cv2
            !! Radial traction is prescribed
            S_U(smtot  , sph_diag-5, j) = -1.0_dp / sqrt(3.0_dp)
            S_U(smtot  , sph_diag-4, j) = sc_b%Bound_sigmarad_jm2
            S_U(smtot  , sph_diag-3, j) = sc_b%Bound_sigmarad_j2
            S_U(smtot  , sph_diag-2, j) = sc_b%Bound_sigmarad_jp2

            T_U(tmtot, tor_diag-3, j) = 0.5_dp
            T_U(tmtot, tor_diag-0, j) = 0.5_dp
        elseif (SWITCH_top_boundary_cond == 'fsv') then
            !! Tangnetial traction zero
            S_U(smtot-1, sph_diag-3, j) = sc_b%Bound_sigmatan_jm2
            S_U(smtot-1, sph_diag-2, j) = sc_b%Bound_sigmatan_j2
            S_U(smtot-1, sph_diag-1, j) = sc_b%Bound_sigmatan_jp2
            !! Radial velocity prescribed
            cv1 =  sqrt (jj         / (2.0_dp*jj + 1.0_dp))
            cv2 = -sqrt((jj+1.0_dp) / (2.0_dp*jj + 1.0_dp))
            S_U(smtot  , sph_diag-8, j) = 0.5_dp * cv1
            S_U(smtot  , sph_diag-7, j) = 0.5_dp * cv2
            S_U(smtot  , sph_diag-1, j) = 0.5_dp * cv1
            S_U(smtot  , sph_diag+0, j) = 0.5_dp * cv2

            T_U(tmtot, tor_diag-2, j) = tc_b%Bound_j_s_jm1
            T_U(tmtot, tor_diag-1, j) = tc_b%Bound_j_s_jp1
        else
            print*, 'Invalid boundary condition at the top boundary...'
            stop
        end if
    end if
end subroutine top_bottom_boundaries



end module mod_matrix