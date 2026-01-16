module mod_solver
    use mod_variables
    use mod_rhs
    use mod_matrix
    use mod_SHTns
    implicit none
    
contains

subroutine solver(t2)
    implicit none
    integer  :: j, jmind
    real(dp) :: t2
    complex(dp), allocatable :: sph(:), tor(:)

    !$OMP PARALLEL PRIVATE(j) PRIVATE(sph, tor)
    allocate(sph(size(RHS%spheriodal)), tor(size(RHS%torodial)))
    !$OMP DO 
    do jmind = 2, compute_scalar_harm
        j = jmindx(jmind,2);
        call generate_rhs(t2, jmind, sph(:), tor(:))
        call lu_solve(S_U(:,:,j), smtot, sph_lo_diag, sph_hi_diag, smtot, sph_band, S_L(:,:,j), sph_lo_diag, S_I(:,j), sph(:))
        call lu_solve(T_U(:,:,j), tmtot, tor_lo_diag, tor_hi_diag, tmtot, tor_band, T_L(:,:,j), tor_lo_diag, T_I(:,j), tor(:))
        call separate_variables(jmind, sph(:), tor(:))
    end do
    !$OMP END DO
    deallocate(sph, tor)
    !$OMP END PARALLEL
end subroutine solver

subroutine separate_variables(jmind, sph, tor)
    implicit none
    integer :: il, ir, jm1, jp1, j00, jms, j, m, jmind
    integer :: nn, smm, tmm
    complex(dp) :: sph(*), tor(*)

    j = jmindx(jmind,2); m = jmindx(jmind,3)
    jm1 = jml(j, m, j-1)
    j00 = jml(j, m, j  )
    jp1 = jml(j, m, j+1)
    jms =  jm(j, m)
    
! !==========VELOCITY & STRESS==========
    smm = 0; tmm = 0
    do il = 1, numoflayers
        nn = layer(il)%n
        do ir = 1, nn
            layer(il)%vel(jm1,ir) = sph(smm + 7*(ir-1) + 1)
            layer(il)%vel(jp1,ir) = sph(smm + 7*(ir-1) + 2)
            layer(il)%vel(j00,ir) = tor(tmm + 3*(ir-1) + 1)
            
            layer(il)%rhoe(ir, jms) = sph(smm + 7*(ir-1)+3)

            layer(il)%str(6*(jms-1)+1,ir) = sph(smm + 7*(ir-1) + 4)
            layer(il)%str(6*(jms-1)+2,ir) = sph(smm + 7*(ir-1) + 5)
            layer(il)%str(6*(jms-1)+4,ir) = sph(smm + 7*(ir-1) + 6)
            layer(il)%str(6*(jms-1)+6,ir) = sph(smm + 7*(ir-1) + 7)

            layer(il)%str(6*(jms-1)+3,ir) = tor(tmm + 3*(ir-1) + 2)
            layer(il)%str(6*(jms-1)+5,ir) = tor(tmm + 3*(ir-1) + 3)
        end do
        smm = sum(layer(1:il)%sm); tmm = sum(layer(1:il)%tm)
        layer(il)%vel(jm1, nn+1) = sph( smm-1)
        layer(il)%vel(jp1, nn+1) = sph( smm  )
        layer(il)%vel(j00, nn+1) = tor( tmm)
    end do
end subroutine separate_variables
    
end module mod_solver