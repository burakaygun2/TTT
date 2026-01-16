module mod_variables
    use mod_parameters
    implicit none

    ! Layer information
    integer :: numoflayers
    integer, allocatable :: layer_ocean(:), layer_ocean_id(:)

    ! Switches
    integer :: SWITCH_self_grav
    integer :: SWITCH_initial_run
    integer :: SWITCH_inertia
    integer :: SWITCH_coriolis
    integer :: SWITCH_vgradv
    integer :: SWITCH_compressible
    integer :: SWITCH_velocity
    integer :: SWITCH_mmtides
    integer :: SWITCH_override_viscosity
    integer :: SWITCH_hyperviscosity
    character(1) :: SWITCH_solid_layer      ! Elastic: e; Maxwell: m; Andrade: a
    character(1) :: SWITCH_ocean_status     ! Viscous: v; Hydrostatic: h
    character(1) :: SWITCH_hydrostatic_core ! Hydrostatic: h; Nonhydrostatic: n
    character(1) :: SWITCH_one_layer        ! Single layer: s; Multi layer: m
    character(1) :: SWITCH_nondim           ! Nondim: n; Dim: d
    character(len=3) :: SWITCH_bottom_boundary_cond    ! No-slip: ns; free-slip: fs; boundary forcing
    character(len=4) :: SWITCH_top_boundary_cond    ! No-slip: ns; free-slip: fs; boundary forcing

    


    ! Spherical harmonic degrees
    integer :: jstart, jfinal, jmax, old_jmax, max_jmax, init_jmax, mmax, jmax_potential
    integer :: total_scalar_harm, compute_scalar_harm
    integer :: total_vector_harm, compute_vector_harm

    ! Time stepping
    integer  :: num_time_steps, pc, signal
    real(dp) :: time_step, dtmin_r, dtmin_t, change_time_step, t, min_dr, weight_sg, CN_theta
    real(dp) :: peaks_r22(1:3), peak_time_r22(1:3)
    integer  :: consecutive_count
    real(dp) :: bottom_grav_acc          ! Bottom gravitational acceleration for single layer model
    character(len=6) :: bf_location             ! top or bottom boundary forcing location

    ! Parameters to be read from file
    real(dp) :: ecc, obl, ang_vel, rotation_period, final_period, moi, Ocean_viscosity, orb
    character(1) :: unit_period

    ! Radial profiles from files
    real(dp), allocatable :: temp_radius(:), temp_density(:), temp_shear_m(:), temp_bulk_m(:), temp_visc(:)
    integer :: N_initial, num_oceans

    ! Auxiliary
    integer :: integral_on
    logical :: constant_density = .FALSE.
    real(dp), allocatable :: rad(:)

    ! Layer type definition
    type :: spherical_layers
        integer :: n, sm, tm, lbound, ubound
        real(dp) :: eta, rho
        real(dp), allocatable :: r_i(:), r_l(:), rho_i(:), rho_l(:), rho_r_i(:), rho_r_l(:)
        real(dp), allocatable :: g_i(:), g_l(:), mass_i(:), mass_l(:), mu_i(:), eta_i(:), K_i(:)
        complex(dp), allocatable :: rhoe(:,:), str(:,:), vel(:,:), pvel(:,:), pvelt(:,:)
        complex(dp), allocatable :: dis(:,:), ur_i(:,:), ut_i(:,:), pot(:,:), grad_pot(:,:), disl(:,:)
        complex(dp), allocatable :: fc(:,:), fc0(:,:), vgradv(:,:), vgradv0(:,:)
        complex(dp), allocatable :: maxwell_term(:,:), andrade_term(:,:), str_a(:,:,:)
        complex(dp), allocatable :: lap_density(:,:), rhoe_l(:,:), divergence_u20(:)
        complex(dp), allocatable :: div_tensor(:,:), div_v(:,:), grad_v(:,:), Qjm(:,:), Sjm(:,:), Tjm(:,:)
        complex(dp), allocatable :: dSjm(:,:), dTjm(:,:)
        complex(dp), allocatable :: curlv_Qjm(:,:), curlv_Sjm(:,:), curlv_Tjm(:,:)
        real(dp), allocatable :: andrade_CC(:), spectrum(:), temp_array2(:), spectrum_v(:), dissipation_time(:)
        real(dp), allocatable :: Vr(:,:,:), Vt(:,:,:), Vp(:,:,:), VV(:,:,:)
        real(dp), allocatable :: curlv_Vr(:,:,:), curlv_Vt(:,:,:), curlv_Vp(:,:,:)
        real(dp), allocatable :: cp_curlVr(:,:,:), cp_curlVt(:,:,:), cp_curlVp(:,:,:)
        real(dp) :: dissipation, average_dissip, v_norm
    end type

    type(spherical_layers), allocatable :: layer(:)
    type(spherical_layers) :: lsingle
    integer :: Nlay, Nrad

    ! Matrix variables
    integer :: smtot, tmtot
    integer, allocatable :: jmindx(:,:), jmlindx(:,:), S_I(:,:), T_I(:,:)
    real(dp), allocatable :: S_U(:,:,:), T_U(:,:,:), S_L(:,:,:), T_L(:,:,:)
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: time_ecc(:), ecc_arr(:)
    complex(dp), allocatable :: indpot(:,:), topur(:,:)


    ! RHS array type
    type :: rhs_array
        complex(dp), allocatable :: spheriodal(:)
        complex(dp), allocatable :: torodial(:)
    end type

    type(rhs_array) :: RHS

    ! Input section
    character(200) :: file_radial_profile, dir_output, body_name
    real(dp) :: alpha, zeta

    ! Love numbers
    complex(dp), allocatable :: love_k(:), forcing_k(:), love_h(:), love_l(:)
    real(dp) :: maxk20, maxrk22, maxik22, maxrk21, maxik21
    real(dp) :: maxh20, maxrh22, maxih22, maxrh21, maxih21
    real(dp) :: maxl20, maxrl22, maxil22, maxrl21, maxil21
    real(dp) :: k20, real_k22, imag_k22, h20, real_h22, imag_h22, l20, real_l22, imag_l22
    real(dp) :: gp20_0, gp20_1, gp20_2, rgp22_0, rgp22_1, rgp22_2, igp22_0, igp22_1, igp22_2
    real(dp) :: peak_t20_1, peak_t20_2, peak_rt22_1, peak_rt22_2, peak_it22_1, peak_it22_2, pt_d_1, pt_d_2
    real(dp) :: dissipation_0, dissipation_1, dissipation_2, avrg_dissip
    real(dp), allocatable :: kh(:), khp(:,:), kh_per(:,:)
    real(dp), allocatable :: gp20(:), rgp22(:), igp22(:)
    real(dp), allocatable :: ur20(:), rur22(:), iur22(:)
    real(dp), allocatable :: ut20(:), rut22(:), iut22(:)
    real(dp) :: k20_p(1:Max_num_period), rk22_p(1:Max_num_period), ik22_p(1:Max_num_period)
    real(dp) :: h20_p(1:Max_num_period), rh22_p(1:Max_num_period), ih22_p(1:Max_num_period)
    real(dp) :: l20_p(1:Max_num_period), rl22_p(1:Max_num_period), il22_p(1:Max_num_period)
    real(dp) :: dissip_p(1:Max_num_period)

    ! File handling
    integer :: iu(1:20), status
    integer :: WRITE_density, WRITE_velocity, WRITE_stress, WRITE_displacement
    integer :: checkpoint, start_output, write_output, start_simulation, start_cp

    ! Grid and spatial variables
    type(spherical_layers), allocatable :: layer_grid(:)
    real(dp), allocatable :: theta(:), phi(:), rho_phi(:,:)
    complex(dp), allocatable :: velocity_ocean_m2(:,:)

contains

    pure integer function jml2(j, m, l)
        implicit none
        integer, intent(in) :: j, m, l
        if (j == 0) then
            jml2 = 1
        else if (j == 1) then
            jml2 = l + 2
        else if (m == 0) then
            jml2 = (2*j - 2)*3 + l - j
        else
            jml2 = (2*j - 1)*3 + l - j
        end if
    end function jml2

    pure integer function jm2(j, m)
        implicit none
        integer, intent(in) :: j, m
        if (j == 0) then
            jm2 = 1
        else if (j == 1) then
            jm2 = 2
        else if (m == 0) then
            jm2 = 2*j - 1
        else
            jm2 = 2*j
        end if
    end function jm2

    pure integer function dirac(a, b)
        implicit none
        integer, intent(in) :: a, b
        dirac = merge(1, 0, a == b)
    end function dirac

    pure integer function jm(j, m)
        implicit none
        integer, intent(in) :: j, m
        if (j <= mmax) then
            jm = j*(j+1)/2 + m + 1
        else
            jm = (mmax+1)*j - ((mmax-1)*mmax/2 - 1) + m - mmax
        end if
    end function jm

    pure integer function jml(j, m, l)
        implicit none
        integer, intent(in) :: j, m, l
        if (j <= mmax) then
            jml = 3 * (j*(j+1)/2 + m) + l - j
        else
            jml = 3*((mmax+1)*j - ((mmax-1)*mmax/2 - 1) + m - mmax) + l - j - 3
        end if
    end function jml

end module mod_variables