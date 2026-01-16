module mod_matrix_coeffs
    use mod_parameters
    use mod_variables
    implicit none
    type :: sph_coeffs
        real(dp) :: mass_bal_rhojm, mass_bal_vjm1, mass_bal_vjp1
        real(dp) :: const_p_vjm1(1:2), const_p_vjp1(1:2), const_p_sj0
        real(dp) :: Mom__jm1_vjm1, Mom__jm1_sj0(1:2), Mom__jm1_sjm2(1:2), Mom__jm1_sj2(1:2), Mom__jm1_rho
        real(dp) :: Mom__jp1_vjp1, Mom__jp1_sj0(1:2), Mom__jp1_sjp2(1:2), Mom__jp1_sj2(1:2), Mom__jp1_rho
        real(dp) :: Const_jm2_vjm1(1:2)
        real(dp) :: Const_j2__vjm1(1:2), Const_j2__vjp1(1:2)
        real(dp) :: Const_jp2_vjp1(1:2)
    end type sph_coeffs

    type :: sph_boundary_coeffs
        real(dp) :: Bound_jm1_s_0, Bound_jm1_s_j2, Bound_jm1_s_jm2
        real(dp) :: Bound_jm1_vjm1, Bound_jm1_vjp1
        real(dp) :: Bound_jp1_s_0, Bound_jp1_s_j2, Bound_jp1_s_jp2
        real(dp) :: Bound_jp1_vjm1, Bound_jp1_vjp1
        real(dp) ::                    Bound_sigmatan_j2, Bound_sigmatan_jm2, Bound_sigmatan_jp2
        real(dp) :: Bound_sigmarad_j0, Bound_sigmarad_j2, Bound_sigmarad_jm2, Bound_sigmarad_jp2
    end type sph_boundary_coeffs

    type :: tor_coeffs
        real(dp) :: Mom__j_v, Mom__j_s_jm1(1:2), Mom__j_s_jp1(1:2)
        real(dp) :: Const_jm1_vj(1:2), Const_jp1_vj(1:2)
    end type tor_coeffs

    type :: tor_boundary_coeffs
        real(dp) :: Bound_j_s_jm1, Bound_j_s_jp1, Bound_j_vj
    end type tor_boundary_coeffs

    type(sph_coeffs)   :: sc
    type(sph_boundary_coeffs) :: sc_b
    type(tor_coeffs)     :: tc
    type(tor_boundary_coeffs) :: tc_b
contains
    

subroutine matrix_coeffs(j)
    implicit none
    real(dp) :: j
    
    call spheriodal_coeffs(j)
    call spheriodal_boundary_coeffs(j)
    call toroidal_coeffs(j)
    call toroidal_boundary_coefss(j)
end subroutine matrix_coeffs

subroutine spheriodal_coeffs(j)
    implicit none
    real(dp) :: j
    !! Balance of mass rho_t = -div(rho0 v)
    !! Below the coefficients of rho_t and vr
    sc%mass_bal_rhojm = 1.0_dp / time_step
    sc%mass_bal_vjm1  =  sqrt(  j        / (2.0_dp * j + 1.0_dp))
    sc%mass_bal_vjp1  = -sqrt( (j+1.0_dp) / (2.0_dp * j + 1.0_dp))

    !! Constitutive - pressure if solid   p_t - vr rho0 g0 = -K div(v)
    !!                         if liquid  p   = -K * rho/rho0 + zeta * div(v)
    !!                            last term in the liquid case may be negligible
    !! below coefficients are just give for div(v)
    sc%const_p_vjm1(1) =  sqrt(  j        / (2.0_dp * j + 1.0_dp))
    sc%const_p_vjp1(1) = -sqrt( (j+1.0_dp) / (2.0_dp * j + 1.0_dp))
    sc%const_p_vjm1(2) = -(j-1.0_dp)*sc%const_p_vjm1(1)
    sc%const_p_vjp1(2) =  (j+2.0_dp)*sc%const_p_vjp1(1)
    !! Balance of momentum rho0 v_t = div(sigma) - rho * g *er - rho0(nabla_V) - 2rho0 omega x v
    !!                         |         implicit            |  |          explicit/RHS         |
    !! below just the coefficients of implicit part
    ! !Y_{JM}^{J-1}
    sc%Mom__jm1_sj0(1)  = -sqrt( j / (3.0_dp * (2.0_dp * j + 1.0_dp)))
    sc%Mom__jm1_sjm2(1) =  sqrt((j - 1.0_dp) / (2.0_dp * j - 1.0_dp))
    sc%Mom__jm1_sj2(1)  = -sqrt( (j + 1.0_dp)*(2.0_dp * j + 3.0_dp) /(6._dp*(2._dp*j - 1._dp)*(2._dp*j + 1._dp) ))
    sc%Mom__jm1_sj0(2)  =  sc%Mom__jm1_sj0(1)  * (j + 1.0_dp)
    sc%Mom__jm1_sjm2(2) = -sc%Mom__jm1_sjm2(1) * (j - 2.0_dp)
    sc%Mom__jm1_sj2(2)  =  sc%Mom__jm1_sj2(1)  * (j + 1.0_dp)
    sc%Mom__jm1_rho     = -sqrt( j       /(2.0_dp * j + 1.0_dp))

    ! !Y_{JM}^{J+1}
    sc%Mom__jp1_sj0(1)   =  sqrt((j + 1.0_dp)/((3._dp*(2._dp*j + 1._dp))))
    sc%Mom__jp1_sjp2(1)  = -sqrt((j +2._dp) / (2._dp*j + 3._dp))
    sc%Mom__jp1_sj2(1)   =  sqrt(( j*(2._dp*j - 1._dp) )/(6._dp*(2._dp*j + 3._dp)*(2._dp*j + 1._dp)))
    sc%Mom__jp1_sj0(2)   = -sc%Mom__jp1_sj0(1)  * j
    sc%Mom__jp1_sjp2(2)  =  sc%Mom__jp1_sjp2(1) * (j + 3.0_dp)
    sc%Mom__jp1_sj2(2)   = -sc%Mom__jp1_sj2(1)  * j
    sc%Mom__jp1_rho      = -(-sqrt((j+1.0_dp)/(2.0_dp * j + 1.0_dp)))

    !!Constitutive - Deviatoric (sigma^d) / eta + (sigma^d)_t / mu  = (nabla_v + nabla_v^T)
    ! ! Y_{JM}^{J-2,2}
    sc%Const_jm2_vjm1(1) = -sqrt((j - 1.0_dp) / (2.0_dp * j - 1.0_dp))
    sc%Const_jm2_vjm1(2) =  sc%Const_jm2_vjm1(1) * j

    ! ! Y_{JM}^{J  ,2}
    sc%Const_j2__vjm1(1) =  sqrt( ((j + 1._dp)*(2._dp*j + 3._dp))/(6._dp*(2._dp*j - 1._dp)*(2._dp*j + 1._dp)))
    sc%Const_j2__vjp1(1) = -sqrt((j*(2._dp*j - 1._dp)) / (6._dp*(2._dp*j + 3._dp)*(2._dp*j + 1._dp)))
    sc%Const_j2__vjm1(2) = -sc%Const_j2__vjm1(1) * (j - 1.0_dp)
    sc%Const_j2__vjp1(2) =  sc%Const_j2__vjp1(1) * (j + 2.0_dp)

    ! ! Y_{JM}^{J+2,2}
    sc%Const_jp2_vjp1(1) =  sqrt( (j + 2.0_dp)  / (2.0_dp * j + 3.0_dp) )
    sc%Const_jp2_vjp1(2) = -(j + 1.0_dp) * sc%Const_jp2_vjp1(1)
end subroutine spheriodal_coeffs

subroutine spheriodal_boundary_coeffs(j)
    implicit none
    real(dp) :: j
    ! !BOUNDARY CONDITIONS Y_{JM}^{J-1}
    sc_b%Bound_jm1_s_0   = -sqrt(j / (3._dp*(2._dp*j +1._dp)) )
    sc_b%Bound_jm1_s_j2  = -sqrt( ((j + 1._dp)*(2._dp*j + 3._dp)) / (6._dp*(2._dp*j + 1._dp)*(2._dp*j -1.0_dp)) )
    sc_b%Bound_jm1_s_jm2 =  sqrt( (j - 1._dp)/(2._dp*j -1._dp) )

    sc_b%Bound_jm1_vjm1 =                        (j / (2.0_dp * j + 1.0_dp))
    sc_b%Bound_jm1_vjp1 = - sqrt(j * (j + 1.0_dp) ) / (2.0_dp * j + 1.0_dp)
    ! !BOUNDARY CONDITIONS Y_{JM}^{J+1}
    sc_b%Bound_jp1_s_0   =  sqrt( (j + 1._dp) / (3._dp*(2._dp*j + 1._dp)) )
    sc_b%Bound_jp1_s_j2  =  sqrt( (  j*(2._dp*j - 1._dp)) / (6._dp*(2._dp*j + 1._dp)*(2._dp*j + 3._dp)) )
    sc_b%Bound_jp1_s_jp2 = -sqrt( (  j + 2._dp) / (2._dp*j + 3._dp) )

    sc_b%Bound_jp1_vjm1 = - sqrt(j * (j + 1.0_dp) ) / (2.0_dp * j + 1.0_dp)
    sc_b%Bound_jp1_vjp1 =             (j + 1.0_dp)   / (2.0_dp * j + 1.0_dp)

    sc_b%Bound_sigmatan_jm2 =  (j+1.0_dp)*sqrt((j-1.0_dp)/(2.0_dp*j-1.0_dp))/(2.0_dp*j+1.0_dp) 
    sc_b%Bound_sigmatan_j2  = -sqrt(3*(j+1.0_dp)/(2.0_dp*(2.0_dp*j-1.0_dp)*(2.0_dp*j+1.0_dp)*(2.0_dp*j+3.0_dp)))
    sc_b%Bound_sigmatan_jp2 = -sqrt(j*(j+1.0_dp)*(j+2.0_dp)/(2.0_dp*j+3.0_dp))/(2.0_dp*j+1.0_dp)

    sc_b%Bound_sigmarad_j2  = -2.0_dp*sqrt(j*(j+1.0_dp))/sqrt(6*(2.0_dp*j-1.0_dp)*(2.0_dp*j+3.0_dp))
    sc_b%Bound_sigmarad_jm2 =  sqrt(         j*(j-1.0_dp)/((2.0_dp*j-1.0_dp)*(2.0_dp*j+1.0_dp)))
    sc_b%Bound_sigmarad_jp2 =  sqrt((j+1.0_dp)*(j+2.0_dp)/((2.0_dp*j+1.0_dp)*(2.0_dp*j+3.0_dp)))
end subroutine spheriodal_boundary_coeffs

subroutine toroidal_coeffs(j)
    implicit none
    real(dp) :: j

    tc%Mom__j_s_jm1(1) =  sqrt((j - 1.0_dp)/(2.0_dp*(2.0_dp * j + 1.0_dp)))
    tc%Mom__j_s_jm1(2) = -(j - 1.0_dp) * tc%Mom__j_s_jm1(1)
    tc%Mom__j_s_jp1(1) = -sqrt((j + 2.0_dp)/(2.0_dp*(2.0_dp * j + 1.0_dp)))
    tc%Mom__j_s_jp1(2) =  (j + 2.0_dp) * tc%Mom__j_s_jp1(1)

    tc%Const_jm1_vj(1) = -sqrt((j - 1.0_dp)/(2.0_dp*(2.0_dp * j + 1.0_dp)))
    tc%Const_jm1_vj(2) =  (j + 1.0_dp) * tc%Const_jm1_vj(1)

    tc%Const_jp1_vj(1) =  sqrt((j + 2.0_dp)/(2.0_dp*(2.0_dp * j + 1.0_dp)))
    tc%Const_jp1_vj(2) = - j * tc%Const_jp1_vj(1)
end subroutine toroidal_coeffs

subroutine toroidal_boundary_coefss(j)
    implicit none
    real(dp) :: j
    
    tc_b%Bound_j_s_jm1 =  sqrt((j-1._dp)/(2._dp*(2._dp*j + 1._dp)))
    tc_b%Bound_j_s_jp1 = -sqrt((j+2._dp)/(2._dp*(2._dp*j + 1._dp)))
end subroutine toroidal_boundary_coefss
end module mod_matrix_coeffs