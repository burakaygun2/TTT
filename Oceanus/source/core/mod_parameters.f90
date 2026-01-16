module mod_parameters
    use iso_c_binding
    use iso_fortran_env, only: real64, int64
    implicit none
    !!Constants
    integer    , parameter  :: dp             = real64
    real(dp)   , parameter  :: Pi             = acos(-1.0_dp)
    real(dp)   , parameter  :: Grav_Const     = 6.67408_dp * (1e-11_dp)
    !! Conversions
    real(dp)   , parameter  :: gcm3_to_kgm3   = 1000.0_dp
    real(dp)   , parameter  :: kgm3_to_gcm3   = 1e-3_dp
    real(dp)   , parameter  :: m_to_km        = 1e-3_dp
    real(dp)   , parameter  :: km_to_m        = 1e3_dp
    real(dp)   , parameter  :: GPa_to_Pa      = 1e9_dp
    real(dp)   , parameter  :: Pa_to_GPa      = 1e-9_dp
    real(dp)   , parameter  :: W_to_GW        = 1e-9_dp
    real(dp)   , parameter  :: W_to_TW        = 1e-12_dp
    real(dp)   , parameter  :: rad_to_deg     = 180.0_dp/Pi
    real(dp)   , parameter  :: deg_to_rad     = Pi/180.0_dp
    real(dp)   , parameter  :: pjmstr         = 1.0_dp/(3.0_dp)**0.5
    real(dp)   ,  parameter :: day2sec        = 86400.0_dp
    !!Zeros
    complex(dp), parameter  :: Im             = complex(0.0_dp, 1.0_dp)
    complex(dp), parameter  :: dzero          = complex(0.0_dp, 0.0_dp)
    real(dp)   , parameter  :: zero           = 0.0_dp
    !! Numerical parameters
    integer   , parameter  :: Max_num_period  = 2000
    integer   , parameter  :: sph_lo_diag     = 12
    integer   , parameter  :: sph_hi_diag     = 9
    integer   , parameter  :: tor_lo_diag     = 3
    integer   , parameter  :: tor_hi_diag     = 4
    integer   , parameter  :: sph_band        = sph_lo_diag + sph_hi_diag + 1
    integer   , parameter  :: tor_band        = tor_lo_diag + tor_hi_diag + 1
    integer   , parameter  :: sph_diag        = sph_lo_diag + 1
    integer   , parameter  :: tor_diag        = tor_lo_diag + 1
    real(dp)  , parameter  :: CFL             = 1.0_dp

    !! EXPERIMENTAL: Bulk viscosity of the ocean for compressible models
    real(dp)  , parameter  :: bulk_vis = 1e-3_dp
    real(dp)  , parameter  :: semimajor_a = 1070400e3_dp
end module mod_parameters