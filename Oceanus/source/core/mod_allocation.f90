module mod_allocation
    use mod_parameters
    use mod_variables
    use mod_plib
    implicit none
contains
    
subroutine alloc_NM
    implicit none

    call delloc_real3(S_U)
    call delloc_real3(S_L)
    call delloc_real3(T_U)
    call delloc_real3(T_L)
    call delloc_int2(S_I)
    call delloc_int2(T_I)

    call alloc_real3(S_U, smtot, sph_band   , Max_jmax)
    call alloc_real3(S_L, smtot, sph_lo_diag, Max_jmax)
    call alloc_real3(T_U, tmtot, tor_band   , Max_jmax)
    call alloc_real3(T_L, tmtot, tor_lo_diag, Max_jmax)
    call alloc_int2(S_I, smtot, Max_jmax)
    call alloc_int2(T_I, tmtot, Max_jmax)
end subroutine alloc_NM

subroutine alloc_RHS
    implicit none
    call alloc_complex1(RHS%spheriodal, smtot)
    call alloc_complex1(RHS%torodial  , tmtot)
end subroutine alloc_RHS


subroutine alloc_vars
    implicit none
    integer :: il, nn

    do il = 1, numoflayers
        nn = layer(il)%n
        call alloc_complex2(layer(il)%rhoe    ,   nn, total_scalar_harm  )
        call alloc_complex2(layer(il)%vel     ,   total_vector_harm, nn+1)
        call alloc_complex2(layer(il)%str     , 6*total_scalar_harm, nn  )
        call alloc_complex2(layer(il)%pot     ,   total_scalar_harm, nn  )
        call alloc_complex2(layer(il)%grad_pot,   total_vector_harm, nn-1)
        call alloc_complex2(layer(il)%fc      ,   total_vector_harm, nn+1)
        call alloc_complex2(layer(il)%fc0     ,   total_vector_harm, nn+1)
        call alloc_complex2(layer(il)%ur_i    ,   total_scalar_harm, 2)
        call alloc_complex2(layer(il)%ut_i    ,   total_scalar_harm, 2)
        call alloc_complex2(layer(il)%pvel    ,   total_scalar_harm, 2)
        call alloc_complex2(layer(il)%pvelt   ,   total_scalar_harm, 2)
        call alloc_complex2(layer(il)%vgradv  ,   nn-1, total_vector_harm)
        call alloc_complex2(layer(il)%vgradv0 ,   nn-1, total_vector_harm)
        
        call alloc_complex2(layer(il)%div_tensor, nn-1, total_vector_harm)
        call alloc_complex2(layer(il)%grad_v    , nn  ,6*total_scalar_harm)
        call alloc_real1   (layer(il)%spectrum  , Max_jmax)
        call alloc_real1   (layer(il)%spectrum_v, Max_jmax)
        call alloc_complex2(layer(il)%Qjm, total_scalar_harm, nn)
        call alloc_complex2(layer(il)%Sjm, total_scalar_harm, nn)
        call alloc_complex2(layer(il)%Tjm, total_scalar_harm, nn)

        call alloc_real1   (layer(il)%dissipation_time, num_time_steps)
    end do

    call alloc_complex1(love_k   , total_scalar_harm)
    call alloc_complex1(love_h   , total_scalar_harm)
    call alloc_complex1(love_l   , total_scalar_harm)
    call alloc_complex1(forcing_k, total_scalar_harm)
    call alloc_real1   (time     , num_time_steps)
    call alloc_complex2(indpot   , jm(2,2), num_time_steps)
    call alloc_complex2(topur    , jm(2,2), num_time_steps)
    
    ! call alloc_real1(kh,3*10)
    ! call alloc_real2(khp, 3*10, num_time_steps)
    ! call alloc_real2(kh_per, 3*10, Max_num_period)
    ! call alloc_real1(gp20 , num_time_steps)
    ! call alloc_real1(rgp22, num_time_steps)
    ! call alloc_real1(igp22, num_time_steps)
    ! call alloc_real1(ur20 , num_time_steps)
    ! call alloc_real1(rur22, num_time_steps)
    ! call alloc_real1(iur22, num_time_steps)
    ! call alloc_real1(ut20 , num_time_steps)
    ! call alloc_real1(rut22, num_time_steps)
    ! call alloc_real1(iut22, num_time_steps)
end subroutine alloc_vars

subroutine alloc_rads
    implicit none
    integer :: il, nn

    do il = 1, numoflayers
        nn = layer(il)%n
        call delloc_real1(layer(il)%r_i    ); call delloc_real1(layer(il)%r_l    );
        call delloc_real1(layer(il)%rho_i  ); call delloc_real1(layer(il)%rho_l  );
        call delloc_real1(layer(il)%rho_r_i); call delloc_real1(layer(il)%rho_r_l);
        call delloc_real1(layer(il)%g_i    ); call delloc_real1(layer(il)%g_l    );
        call delloc_real1(layer(il)%mass_i ); call delloc_real1(layer(il)%mass_l );
        call delloc_real1(layer(il)%mu_i   );
        call delloc_real1(layer(il)%eta_i  );
        call delloc_real1(layer(il)%K_i    );

        call alloc_real1(layer(il)%r_i    , nn); call alloc_real1(layer(il)%r_l    , nn+1);
        call alloc_real1(layer(il)%rho_i  , nn); call alloc_real1(layer(il)%rho_l  , nn+1);
        call alloc_real1(layer(il)%rho_r_i, nn); call alloc_real1(layer(il)%rho_r_l, nn+1);
        call alloc_real1(layer(il)%g_i    , nn); call alloc_real1(layer(il)%g_l    , nn+1);
        call alloc_real1(layer(il)%mass_i , nn); call alloc_real1(layer(il)%mass_l , nn+1);
        call alloc_real1(layer(il)%mu_i   , nn);
        call alloc_real1(layer(il)%eta_i  , nn);
        call alloc_real1(layer(il)%K_i    , nn);
    end do
end subroutine alloc_rads


subroutine alloc_single_layer
    implicit none

    call alloc_real1(lsingle%r_i, lsingle%n  )
    call alloc_real1(lsingle%r_l, lsingle%n+1)
    
end subroutine alloc_single_layer

subroutine alloc_rads_lsingle
    implicit none
    
    call alloc_real1(lsingle%rho_i  , lsingle%n  )
    call alloc_real1(lsingle%rho_l  , lsingle%n+1)
    call alloc_real1(lsingle%g_i    , lsingle%n  )

    call alloc_real1(lsingle%eta_i  , lsingle%n)
end subroutine alloc_rads_lsingle

subroutine alloc_vars_lsingle
    implicit none
    call alloc_complex2(lsingle%rhoe     , lsingle%n, total_scalar_harm) 
    call alloc_complex2(lsingle%vel      , total_vector_harm, lsingle%n+1)
    call alloc_complex2(lsingle%str      , 6*total_scalar_harm, lsingle%n)
    call alloc_complex2(lsingle%pot      , total_scalar_harm, lsingle%n)
    call alloc_complex2(lsingle%grad_pot , total_vector_harm, lsingle%n-1)
    call alloc_complex2(lsingle%fc       , total_vector_harm, lsingle%n+1)
    call alloc_complex2(lsingle%fc0      , total_vector_harm, lsingle%n+1)
    call alloc_complex2(lsingle%ur_i     , total_scalar_harm, 2)
    
end subroutine alloc_vars_lsingle


subroutine alloc_grid
    implicit none
    
    call delloc_real2(rho_phi)
    call delloc_real1(theta)
    call delloc_real1(phi)

    call  alloc_real2(rho_phi, layer(numoflayers)%n-1, 360)
    call  alloc_real1(theta, 360)
    call  alloc_real1(phi,   360)
end subroutine alloc_grid


end module mod_allocation
