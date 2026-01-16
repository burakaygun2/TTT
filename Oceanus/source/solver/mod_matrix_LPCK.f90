module mod_matrix_LPCK
    use mod_allocation
    use mod_variables
    use mod_parameters
    use mod_matrix_coeffs
    use mod_LU
    implicit none
    real(dp), allocatable :: AA(:,:), TAA(:,:)
contains

subroutine create_matrix_LPCK(j)
    implicit none
    integer :: j!, info
    real(dp) :: jj!, d,d2


    jj = real(j, dp)
    call matrix_coeffs(jj)
    call generate_full_matrix
    call convert_to_band(j)
    if(SWITCH_compressible==0)then

        ! call generate_matrix_incompressible(j)
        ! call internal_boundaries(j)
        ! call top_bottom_boundaries(j)
    elseif(SWITCH_compressible==1)then
        ! call generate_matrix_compressible(j)
        ! call internal_boundaries(j)
        ! call top_bottom_boundaries(j)
    else
        stop "For incompressible: 0 or Compressible: 1..."
    end if
    ! call dgbtrf(smtot, smtot, skl, sku, S_LPCK(:,:,j), Sldab, Spiv(:,j), info)
    ! call dgbtrf(tmtot, tmtot, tor_lo_diag, tor_up_diag, T_LPCK, Tldab, Tpiv, info)
end subroutine create_matrix_LPCK


subroutine generate_full_matrix
    implicit none
    integer :: il, ir, smm, tmm, row, nn
    real(dp) :: hl, hup, hlo, hup_r, hlo_r
    real(dp) :: hi, hup2, hlo2, hup2_r, hlo2_r
    real(dp) :: drho, g, rho

    allocate(AA(smtot, smtot), TAA(tmtot,tmtot))
    AA = 0.0_dp; TAA = 0.0_dp
    smm = 0; tmm = 0
    do il = 1, numoflayers
        do ir = 1, layer(il)%n
            hl  =  layer(il)%r_l(ir+1) - layer(il)%r_l(ir)
            hup = (layer(il)%r_l(ir+1) - layer(il)%r_i(ir)) / hl;           hup_r = hup/layer(il)%r_i(ir)
            hlo = (layer(il)%r_i(ir)   - layer(il)%r_l(ir)) / hl;           hlo_r = hlo/layer(il)%r_i(ir)
            row = smm + 7*(ir-1)+1
            AA(row, row+2) = 1.0_dp !!Balance of mass rho_jm ! for incompressible constant density model this is trivial solution

            row = smm + 7*(ir-1)+2
            AA(row, row-1) = -sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hup_r !!div(v) = 0 
            AA(row, row+0) = -sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hup_r !!div(v) = 0 
            AA(row, row+6) =  sc%const_p_vjm1(1)/hl + sc%const_p_vjm1(2)*hlo_r !!div(v) = 0 
            AA(row, row+7) =  sc%const_p_vjp1(1)/hl + sc%const_p_vjp1(2)*hlo_r !!div(v) = 0 
            !! CONSTITUTIVE EQUATION Y_{JM}^{J-2,2} 
            row = smm + 7*(ir-1)+5
            AA(row, row-4) = (-sc%Const_jm2_vjm1(1)/hl + sc%Const_jm2_vjm1(2)*hup_r)
            AA(row, row+0) = 1.0d0/(2.0d0*layer(il)%eta_i(ir)) + 1.0d0 / (2.0d0 * layer(il)%mu_i(ir)*time_step) 
            AA(row, row+3) = ( sc%Const_jm2_vjm1(1)/hl + sc%Const_jm2_vjm1(2)*hlo_r)
            !! CONSTITUTIVE EQUATION Y_{JM}^{J  ,2}
            row = smm + 7*(ir-1)+6
            AA(row, row-5) = (-sc%Const_j2__vjm1(1)/hl + sc%Const_j2__vjm1(2)*hup_r)
            AA(row, row-4) = (-sc%Const_j2__vjp1(1)/hl + sc%Const_j2__vjp1(2)*hup_r)
            AA(row, row+0) = 1.0d0/(2.0d0*layer(il)%eta_i(ir)) + 1.0d0 / (2.0d0 * layer(il)%mu_i(ir)*time_step) 
            AA(row, row+2) = ( sc%Const_j2__vjm1(1)/hl + sc%Const_j2__vjm1(2)*hlo_r)
            AA(row, row+3) = ( sc%Const_j2__vjp1(1)/hl + sc%Const_j2__vjp1(2)*hlo_r)
            !! CONSTITUTIVE EQUATION Y_{JM}^{J+2,2}
            row = smm + 7*(ir-1)+7
            AA(row, row-5) = (-sc%Const_jp2_vjp1(1)/hl + sc%Const_jp2_vjp1(2)*hup_r)
            AA(row, row+0) = 1.0d0/(2.0d0*layer(il)%eta_i(ir)) + 1.0d0 / (2.0d0 * layer(il)%mu_i(ir)*time_step) 
            AA(row, row+2) = ( sc%Const_jp2_vjp1(1)/hl + sc%Const_jp2_vjp1(2)*hlo_r)
            if(SWITCH_ocean_status == 'v')then
                ! !CONSTITUTIVE EQUATION Y_{JM}^{J-1,2} 
                row = tmm+3*(ir-1)+2
                TAA(row, row-1) = (-tc%Const_jm1_vj(1)/hl + tc%Const_jm1_vj(2)*hup_r)
                TAA(row, row+0) = 1.0d0 / (2.0d0 * layer(il)%mu_i(ir)*time_step) + 1.0d0/(2.0d0*layer(il)%eta_i(ir))
                TAA(row, row+2) = ( tc%Const_jm1_vj(1)/hl + tc%Const_jm1_vj(2)*hlo_r)

                ! !CONSTITUTIVE EQUATION Y_{JM}^{J+1,2}
                row = tmm+3*(ir-1)+3
                TAA(row, row-2) = (-tc%Const_jp1_vj(1)/hl + tc%Const_jp1_vj(2)*hup_r)
                TAA(row, row+0) = 1.0d0 / (2.0d0 * layer(il)%mu_i(ir)*time_step) + 1.0d0/(2.0d0*layer(il)%eta_i(ir))
                TAA(row, row+1) = ( tc%Const_jp1_vj(1)/hl + tc%Const_jp1_vj(2)*hlo_r)
            end if
        end do
        do ir = 2, layer(il)%n
            hi   =   layer(il)%r_i(ir)   - layer(il)%r_i(ir-1)
            hup2 =  (layer(il)%r_i(ir)   - layer(il)%r_l(  ir)) / hi;           hup2_r = hup2/layer(il)%r_l(ir)
            hlo2 =  (layer(il)%r_l(ir)   - layer(il)%r_i(ir-1)) / hi;           hlo2_r = hlo2/layer(il)%r_l(ir)
            !!Y_{JM}^{J-1} COMPONENT
            row = smm+7*(ir-1)+3
            AA(row, row-6) = (CN_theta)*(-sc%Mom__jm1_sj0(1) /hi + sc%Mom__jm1_sj0(2)  * hup2_r)   ! !S_{JM}^{J  ,0}
            AA(row, row-5) = (CN_theta)*(-sc%Mom__jm1_sjm2(1)/hi + sc%Mom__jm1_sjm2(2) * hup2_r)   ! !S_{JM}^{J-2,2}
            AA(row, row-4) = (CN_theta)*(-sc%Mom__jm1_sj2(1) /hi + sc%Mom__jm1_sj2(2)  * hup2_r)   ! !S_{JM}^{J  ,2}
            AA(row, row+1) = (CN_theta)*( sc%Mom__jm1_sj0(1) /hi + sc%Mom__jm1_sj0(2)  * hlo2_r)   ! !S_{JM}^{J  ,0}
            AA(row, row+2) = (CN_theta)*( sc%Mom__jm1_sjm2(1)/hi + sc%Mom__jm1_sjm2(2) * hlo2_r)   ! !S_{JM}^{J-2,2}
            AA(row, row+3) = (CN_theta)*( sc%Mom__jm1_sj2(1) /hi + sc%Mom__jm1_sj2(2)  * hlo2_r)   ! !S_{JM}^{J  ,2}
            AA(row, row-2) =        -SWITCH_inertia*layer(il)%rho_l(ir)/time_step

            !!Y_{JM}^{J+1} COMPONENT
            row = smm+7*(ir-1)+4
            AA(row, row-7) = (CN_theta)*(-sc%Mom__jp1_sj0(1) /hi + sc%Mom__jp1_sj0(2)  * hup2_r)  ! !S_{JM}^{J  ,0}
            AA(row, row-5) = (CN_theta)*(-sc%Mom__jp1_sj2(1) /hi + sc%Mom__jp1_sj2(2)  * hup2_r)  ! !S_{JM}^{J  ,2}
            AA(row, row-4) = (CN_theta)*(-sc%Mom__jp1_sjp2(1)/hi + sc%Mom__jp1_sjp2(2) * hup2_r)  ! !S_{JM}^{J+2,2}
            AA(row, row+0) = (CN_theta)*( sc%Mom__jp1_sj0(1) /hi + sc%Mom__jp1_sj0(2)  * hlo2_r)  ! !S_{JM}^{J  ,0}
            AA(row, row+2) = (CN_theta)*( sc%Mom__jp1_sj2(1) /hi + sc%Mom__jp1_sj2(2)  * hlo2_r)  ! !S_{JM}^{J  ,2}
            AA(row, row+3) = (CN_theta)*( sc%Mom__jp1_sjp2(1)/hi + sc%Mom__jp1_sjp2(2) * hlo2_r)  ! !S_{JM}^{J+2,2}
            AA(row, row-2) =        -SWITCH_inertia*layer(il)%rho_l(ir)/time_step

            if(SWITCH_ocean_status == 'v')then
                row = tmm+3*(ir-1)+1
                TAA(row, row-2) = (CN_theta)*(-tc%Mom__j_s_jm1(1) /hi + tc%Mom__j_s_jm1(2)  * hup2_r)
                TAA(row, row-1) = (CN_theta)*(-tc%Mom__j_s_jp1(1) /hi + tc%Mom__j_s_jp1(2)  * hup2_r)
                TAA(row, row+0) = -SWITCH_inertia*layer(il)%rho_l(ir)/time_step
                TAA(row, row+1) = (CN_theta)*( tc%Mom__j_s_jm1(1) /hi + tc%Mom__j_s_jm1(2)  * hlo2_r)
                TAA(row, row+2) = (CN_theta)*( tc%Mom__j_s_jp1(1) /hi + tc%Mom__j_s_jp1(2)  * hlo2_r)
            end if
        end do
        smm = sum(layer(1:il)%sm); tmm=sum(layer(1:il)%tm)
    end do
    
    do il = 2, numoflayers
        nn = layer(il-1)%n; smm = sum(layer(1:il-1)%sm); tmm = sum(layer(1:il-1)%tm)
        drho  = (layer(il-1)%rho_i(layer(il-1)%n) - layer(il)%rho_i(1))
        g     = layer(il-1)%  g_i(layer(il-1)%n)
        hl    =  layer(il-1)%r_l(nn+1) - layer(il-1)%r_l(nn)
        hup   = (layer(il-1)%r_l(nn+1) - layer(il-1)%r_i(nn))/(hl);            hup_r  = hup/layer(il-1)%r_i(nn)
        hlo   = (layer(il-1)%r_i(nn)   - layer(il-1)%r_l(nn))/(hl);            hlo_r  = hlo/layer(il-1)%r_i(nn)
        ! !CONTIUNITY OF TRACTION_{JM}^{J-1}
        row = smm+3
        AA(row, row-11) =  -hup * sc_b%Bound_jm1_vjm1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J-1}, LOWER
        AA(row, row-10) =  -hup * sc_b%Bound_jm1_vjp1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J-1}, LOWER
        AA(row, row-8 ) =  -sc_b%Bound_jm1_s_0   ! !S_{JM}^{J  ,0}, LOWER
        AA(row, row-7 ) =  -sc_b%Bound_jm1_s_jm2 ! !S_{JM}^{J-2,2}, LOWER
        AA(row, row-6 ) =  -sc_b%Bound_jm1_s_j2  ! !S_{JM}^{J  ,2}, LOWER
        AA(row, row-4 ) =  -hlo * sc_b%Bound_jm1_vjm1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J-1}, LOWER
        AA(row, row-3 ) =  -hlo * sc_b%Bound_jm1_vjp1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J-1}, LOWER
        AA(row, row+1 ) =   sc_b%Bound_jm1_s_0   ! !S_{JM}^{J  ,0}, UPPER
        AA(row, row+2 ) =   sc_b%Bound_jm1_s_jm2 ! !S_{JM}^{J-2,2}, UPPER
        AA(row, row+3 ) =   sc_b%Bound_jm1_s_j2  ! !S_{JM}^{J  ,2}, UPPER
        ! ?==========================================================
        ! ?==========================================================
        ! !CONTIUNITY OF TRACTION_{JM}^{J-1}
        row = smm+4
        AA(row, row-12) =  -hup * sc_b%Bound_jp1_vjm1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J-1}, LOWER
        AA(row, row-11) =  -hup * sc_b%Bound_jp1_vjp1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J+1}, LOWER
        AA(row, row-9 ) =  -sc_b%Bound_jp1_s_0   ! !S_{JM}^{J  ,0}, LOWER
        AA(row, row-7 ) =  -sc_b%Bound_jp1_s_j2  ! !S_{JM}^{J-2,2}, LOWER
        AA(row, row-6 ) =  -sc_b%Bound_jp1_s_jp2 ! !S_{JM}^{J  ,2}, LOWER
        AA(row, row-5 ) =  -hlo * sc_b%Bound_jp1_vjm1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J-1}, LOWER
        AA(row, row-4 ) =  -hlo * sc_b%Bound_jp1_vjp1 * drho * g * 0.5d0 * time_step ! !V_{JM}^{J+1}, LOWER
        AA(row, row+0 ) =   sc_b%Bound_jp1_s_0   ! !S_{JM}^{J  ,0}, UPPER
        AA(row, row+2 ) =   sc_b%Bound_jp1_s_j2  ! !S_{JM}^{J  ,2}, UPPER
        AA(row, row+3 ) =   sc_b%Bound_jp1_s_jp2 ! !S_{JM}^{J+2,2}, UPPER
        row = tmm+1
        TAA(row, row-3) = -tc_b%Bound_j_s_jm1 !core
        TAA(row, row-2) = -tc_b%Bound_j_s_jp1 !core
        TAA(row, row+1) =  tc_b%Bound_j_s_jm1!ocean
        TAA(row, row+2) =  tc_b%Bound_j_s_jp1!ocean
    end do
    do il = 1, numoflayers-1
        smm = sum(layer(1:il)%sm); tmm = sum(layer(1:il)%tm)
        ! !CONTIUNITY OF V_{JM}^{J-1}
        row = smm-1
        AA(row, row-7) = -0.5d0 ! !BELOW LAYER
        AA(row, row+0) = -0.5d0 ! !BELOW LAYER
        AA(row, row+2) =  0.5d0 ! !UPPER LAYER
        AA(row, row+9) =  0.5d0 ! !UPPER LAYER

        ! !CONTIUNITY OF V_{JM}^{J+1}
        row = smm
        AA(row, row-7) = -0.5d0 ! !BELOW LAYER
        AA(row, row+0) = -0.5d0 ! !BELOW LAYER
        AA(row, row+2) =  0.5d0 ! !UPPER LAYER
        AA(row, row+9) =  0.5d0 ! !UPPER LAYER
        if(SWITCH_ocean_status=='v')then
            row = tmm
            TAA(row, row-3) = -0.5d0
            TAA(row, row+0) = -0.5d0
            TAA(row, row+1) =  0.5d0
            TAA(row, row+4) =  0.5d0
        end if
    end do

    TAA(1, 1+0) = 1.0d0

    AA(3, 3-2) = 1.0d0
    AA(4, 4-2) = 1.0d0

    il = numoflayers
    nn = layer(il)%n
    rho =  layer(il)%rho_i(nn)
    g   =  layer(il)%g_i(nn)
    hl  =  layer(il)%r_l(nn+1) - layer(il)%r_l(nn)
    hup = (layer(il)%r_l(nn+1) - layer(il)%r_i(nn))/(hl)
    hlo = (layer(il)%r_i(nn)   - layer(il)%r_l(nn))/(hl)

    ! !Free surface: Te_r = - u_r * rho * g = - (u_i + (v_{i+1} + v_{i})*dt/2) * rho * g
    ! !Traction vector j-1
    row = smtot-1
    AA(row, row-7) = sc_b%Bound_jm1_vjm1 * hup * rho * g * 0.5d0 * time_step ! 
    AA(row, row-6) = sc_b%Bound_jm1_vjp1 * hup * rho * g * 0.5d0 * time_step ! 
    AA(row, row-4) = sc_b%Bound_jm1_s_0    ! !S_{JM}^{J  ,0}
    AA(row, row-3) = sc_b%Bound_jm1_s_jm2  ! !S_{JM}^{J-2,2}
    AA(row, row-2) = sc_b%Bound_jm1_s_j2   ! !S_{JM}^{J  ,2}
    AA(row, row+0) = sc_b%Bound_jm1_vjm1* hlo * rho * g * 0.5d0 * time_step !
    AA(row, row+1) = sc_b%Bound_jm1_vjp1* hlo * rho * g * 0.5d0 * time_step !


    ! !Traction vector j+1
    row = smtot
    AA(row  , row-8) = sc_b%Bound_jp1_vjm1* hup * rho * g * 0.5d0 * time_step !
    AA(row  , row-7) = sc_b%Bound_jp1_vjp1* hup * rho * g * 0.5d0 * time_step !
    AA(row  , row-5) = sc_b%Bound_jp1_s_0   ! !S_{JM}^{J  ,0}
    AA(row  , row-3) = sc_b%Bound_jp1_s_j2  ! !S_{JM}^{J  ,2}
    AA(row  , row-2) = sc_b%Bound_jp1_s_jp2 ! !S_{JM}^{J+2,2}
    AA(row  , row-1) = sc_b%Bound_jp1_vjm1 * hlo * rho * g * 0.5d0 * time_step !
    AA(row  , row+0) = sc_b%Bound_jp1_vjp1 * hlo * rho * g * 0.5d0 * time_step !

    row = tmtot
    TAA(row, row-2) = tc_b%Bound_j_s_jm1
    TAA(row, row-1) = tc_b%Bound_j_s_jp1
end subroutine generate_full_matrix

! subroutine convert_to_band(j)
!     implicit none
!     integer :: i, icol, j
!     do icol = 1, smtot
!         do i = max(1, icol - sku), min(smtot, icol + skl)
!             S_LPCK(skl + sku + 1 + i - icol, icol, j) = AA(i, icol)
!         enddo
!     enddo
!     deallocate(AA, TAA)
! end subroutine convert_to_band


end module mod_matrix_LPCK