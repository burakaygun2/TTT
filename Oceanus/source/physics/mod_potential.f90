module mod_potential
    use mod_parameters
    use mod_variables
    use omp_lib
    implicit none
    type :: other_moons
        real(dp) :: mass, ang_vel, semimajor_a, period, rij, nij, C20, C22, S22, cospij, sinpij
        real(dp) :: cos2pij, sin2pij
    end type other_moons
    type(other_moons) :: Europa, Io, Moon
contains

subroutine induced_potential_initial_run(it2)
    implicit none
    integer    :: m, il, ir, jmi, nn, j, it2
    real(dp)    :: cc, jj, cj, ff, ffi
    complex(dp) :: ind_pot_a, ind_pot_b

    if(it2==1) then
        maxk20  = -1.0_dp
        maxrk21 = -1.0_dp
        maxik21 = -1.0_dp
        maxrk22 = -1.0_dp
        maxik22 = -1.0_dp
    end if
    ind_pot_a = dzero
    ind_pot_b = dzero
    ! do il = 1, numoflayers
    !     layer(il)%pot = dzero
    ! end do
    cc = 4.0_dp * Grav_Const * Pi

    ! !$omp parallel do private(m, ir, jmi, jj, cj, ind_pot) shared(cc)

    do il = 1, numoflayers
        do ir = 1,layer(il)%n
            jj = 2.0_dp; cj = 1.0_dp/(2.0_dp*jj+1.0_dp);
            m  = 0; jmi = jm(int(jj), m)
            call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
            layer(il)%pot(jmi, ir) = (ind_pot_a + ind_pot_b) * layer(il)%r_i(ir) * cc * cj

            m  = 1; jmi = jm(int(jj), m)
            call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
            layer(il)%pot(jmi, ir) = (ind_pot_a + ind_pot_b) * layer(il)%r_i(ir) * cc * cj

            m  = 2; jmi = jm(int(jj), m)
            call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
            layer(il)%pot(jmi, ir) = (ind_pot_a + ind_pot_b) * layer(il)%r_i(ir) * cc * cj
        end do
    end do
    ! !$omp end parallel do

    nn = layer(numoflayers)%n
    do j = 2, 2
        jmi = jm(j, 0)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff  = real(forcing_k(jm(2,0)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxk20) then
            maxk20 = real(layer(numoflayers)%pot(jmi,nn))
            love_k(jmi) = cmplx(maxk20/ff, 0.0, dp)
        end if

        ! jmi = jm(j, 1)
        ! indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ! ff = real(forcing_k(jm(2,1)))
        ! if(real(layer(numoflayers)%pot(jmi,nn))>maxrk21) then
        !     maxrk21 = real(layer(numoflayers)%pot(jmi,nn))
        ! end if
        ! love_k(jmi) = cmplx(maxrk21 / ff, 0.0_dp, dp)

        jmi = jm(j, 2)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff = real(forcing_k(jm(2,2)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxrk22) then
            maxrk22 = real(layer(numoflayers)%pot(jmi,nn))
        end if
        ffi = imag(forcing_k(jm(2,2)))
        if(imag(layer(numoflayers)%pot(jmi,nn))>maxik22) then
            maxik22 = imag(layer(numoflayers)%pot(jmi,nn))
        end if
        love_k(jmi) = cmplx(maxrk22 / ff, maxik22 / ffi, dp)
    end do
    if(mod(it2, num_time_steps) == 0) then
        maxk20  = -1.0_dp
        maxrk21 = -1.0_dp
        maxik21 = -1.0_dp
        maxrk22 = -1.0_dp
        maxik22 = -1.0_dp
    end if
end subroutine induced_potential_initial_run

subroutine induced_potential(it2)
    implicit none
    integer    :: m, il, ir, jmi, nn, j, it2
    real(dp)    :: cc, jj, cj, ff, ffi
    complex(dp) :: ind_pot_a, ind_pot_b

    if(it2==1) then
        maxk20  = -1.0_dp
        maxrk21 = -1.0_dp
        maxik21 = -1.0_dp
        maxrk22 = -1.0_dp
        maxik22 = -1.0_dp
    end if
    ind_pot_a = dzero
    ind_pot_b = dzero
    ! do il = 1, numoflayers
    !     layer(il)%pot = dzero
    ! end do
    cc = 4.0_dp * Grav_Const * Pi

    !$omp parallel do private(m, ir, jmi, jj, cj, ind_pot_a, ind_pot_b) shared(cc)

    do il = 1, numoflayers
        do ir = 1,layer(il)%n
            if (SWITCH_coriolis == 0) then
                jj = 2.0_dp; cj = 1.0_dp/(2.0_dp*jj+1.0_dp);
                m  = 0; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 1; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 2; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj
            else
                jj = 2.0_dp; cj = 1.0_dp/(2.0_dp*jj+1.0_dp);
                m  = 0; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 1; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 2; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj
                
                jj = 4.0_dp; cj = 1.0_dp/(2.0_dp*jj+1.0_dp);
                m  = 0; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 1; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 2; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                jj = 6.0_dp; cj = 1.0_dp/(2.0_dp*jj+1.0_dp);
                m  = 0; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 1; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 2; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                jj = 8.0_dp; cj = 1.0_dp/(2.0_dp*jj+1.0_dp);
                m  = 0; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 1; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj

                m  = 2; jmi = jm(int(jj), m)
                call integrator(layer(il)%r_i(ir), int(jj), m, il, ir, ind_pot_a, ind_pot_b)
                layer(il)%pot(jmi, ir) = (ind_pot_a+ind_pot_b) * layer(il)%r_i(ir) * cc * cj
            end if
        end do
    end do
    !$omp end parallel do

    nn = layer(numoflayers)%n
    do j = 2, 2
        jmi = jm(j, 0)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff  = real(forcing_k(jm(2,0)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxk20) then
            maxk20 = real(layer(numoflayers)%pot(jmi,nn))
            love_k(jmi) = cmplx(maxk20/ff, 0.0, dp)
        end if

        jmi = jm(j, 1)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff = real(forcing_k(jm(2,1)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxrk21) then
            maxrk21 = real(layer(numoflayers)%pot(jmi,nn))
        end if
        love_k(jmi) = cmplx(maxrk21 / ff, 0.0_dp, dp)

        jmi = jm(j, 2)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff = real(forcing_k(jm(2,2)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxrk22) then
            maxrk22 = real(layer(numoflayers)%pot(jmi,nn))
        end if
        ffi = imag(forcing_k(jm(2,2)))
        if(imag(layer(numoflayers)%pot(jmi,nn))>maxik22) then
            maxik22 = imag(layer(numoflayers)%pot(jmi,nn))
        end if
        love_k(jmi) = cmplx(maxrk22 / ff, maxik22 / ffi, dp)
    end do
    if(mod(it2, num_time_steps) == 0) then
        maxk20  = -1.0_dp
        maxrk21 = -1.0_dp
        maxik21 = -1.0_dp
        maxrk22 = -1.0_dp
        maxik22 = -1.0_dp
    end if
end subroutine induced_potential

subroutine potential_at_altidude(r, ind_pot)
    implicit none
    integer    :: j, m, il, ir, jmc, nn
    real(dp)    :: dr, cint0, cint1, r, drho, cj, cc
    complex(dp) :: ind_pot
    !! this needs fixing... the indexing is wrong
    cc = 4.0_dp * Grav_Const * Pi
    do j = 2, min(jfinal, 8), 2
        cj = 1.0_dp/(2.0_dp*real(j,dp)+1.0_dp);
        do m = 0, min(j, 2), 2
            jmc = jm2(j, m)
            ind_pot = 0.0_dp
            do il = 3, numoflayers
                do ir = 1, layer(il)%n-1
                    dr = layer(il)%r_i(ir+1) - layer(il)%r_i(ir)
                    cint0 = (layer(il)%r_i(ir)  /r)  ** (j+2)
                    cint1 = (layer(il)%r_i(ir+1)/r)  ** (j+2)
                    ind_pot = ind_pot + 0.5_dp * dr * (cint0 * layer(il)%rhoe(ir, jmc) + cint1 * layer(il)%rhoe(ir+1, jmc))
                end do
            end do

            do il = 1, numoflayers
                nn      = layer(il)%n
                if(il< numoflayers) drho    = layer(il)%rho_i(nn) - layer(il+1)%rho_i(1)
                if(il==numoflayers) drho    = layer(il)%rho_i(nn)
                ind_pot = ind_pot + drho * (layer(il)%ur_i(nn, jmc)) * (layer(il)%r_i(nn)  /r)  ** (j+2)
            end do
            ind_pot = ind_pot * cj
        end do
    end do
    ind_pot = ind_pot * cc * r

end subroutine potential_at_altidude

subroutine self_gravity
    implicit none
    integer :: j, m, il, ir, jm1, jp1, jmi
    real(dp) :: jj, dc1, dc2, c1, c2, hi, hup, hlo

    do il = 1, numoflayers
        do j = 2, min(jmax, 8), 2
            jj  = real(j, dp)
            dc1 =  sqrt( jj       /(2.0_dp*jj+1.0_dp))
            dc2 = -sqrt((jj+1.0_dp)/(2.0_dp*jj+1.0_dp))
            c1  =  (jj+1.0_dp)*dc1
            c2  = -(jj      )*dc2
            do m = 0, mmax
                jm1 = jml(j,m,j-1)
                jp1 = jml(j,m,j+1)
                jmi = jm (j,m)
                do ir = 1, layer(il)%n-1
                    hi  =   layer(il)%r_i(ir+1)  - layer(il)%r_i(ir)
                    hup =  (layer(il)%r_i(ir+1)  - layer(il)%r_l(ir+1)) / hi
                    hlo =  (layer(il)%r_l(ir+1)  - layer(il)%r_i(ir  )) / hi
                    layer(il)%grad_pot(jm1, ir) = dc1 * (layer(il)%pot(jmi, ir+1)     - layer(il)%pot(jmi, ir)    )/hi + &
                                                   c1 * (layer(il)%pot(jmi, ir+1)*hlo + layer(il)%pot(jmi, ir)*hup)/(layer(il)%r_l(ir+1))
                    layer(il)%grad_pot(jp1, ir) = dc2 * (layer(il)%pot(jmi, ir+1)     - layer(il)%pot(jmi, ir)    )/hi + &
                                                   c2 * (layer(il)%pot(jmi, ir+1)*hlo + layer(il)%pot(jmi, ir)*hup)/(layer(il)%r_l(ir+1))
                end do
            end do
        end do
    end do
end subroutine self_gravity

subroutine ind_pot2_coriolis(it2)
    implicit none
    integer  :: j, m, it2, nn, jmi
    real(dp) :: ff, ffi

    jmax_potential = 10
    if(it2==1) then
        maxk20  = -1.0_dp;  maxh20 = -1.0_dp;  maxl20 = -1.0_dp;
        maxrk21 = -1.0_dp; maxrh21 = -1.0_dp; maxrl21 = -1.0_dp;
        maxik21 = -1.0_dp; maxih21 = -1.0_dp; maxil21 = -1.0_dp;
        maxrk22 = -1.0_dp; maxrh22 = -1.0_dp; maxrl22 = -1.0_dp;
        maxik22 = -1.0_dp; maxih22 = -1.0_dp; maxil22 = -1.0_dp;
    end if

    !!$omp parallel do private(j, m)
    do j = 2, jmax_potential, 2
        do m = 0, mmax
            call induced_potential2(j, m)
        end do
    end do
    !!$omp end parallel do

    nn = layer(numoflayers)%n
    !! following parts are the calculations of the Love numbers for degree 2 only
    do j = 2, 2
        jmi = jm(j, 0)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        topur (jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%ur_i(jmi,2)
        ff  = real(forcing_k(jm(2,0)))
        ! print*, real(layer(numoflayers)%pot(jmi,nn)), maxk20
        if(real(layer(numoflayers)%pot(jmi,nn))>maxk20) then
            maxk20 = real(layer(numoflayers)%pot(jmi,nn))
            maxh20 = real(layer(numoflayers)%ur_i(jmi,2))
            maxl20 = real(layer(numoflayers)%ut_i(jmi,2))
            if(ff>0) then
                love_k(jmi) = cmplx(maxk20/ff, 0.0, dp)
                love_h(jmi) = cmplx(maxh20/(ff/layer(numoflayers)%g_i(nn)), 0.0, dp)
                love_l(jmi) = cmplx(maxl20/(ff*sqrt(30.0_dp/5.0_dp)/layer(numoflayers)%g_i(nn)), 0.0, dp)
            end if
        end if

        ! jmi = jm(j, 1)
        ! indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ! ff = real(forcing_k(jm(2,1)))
        ! if(real(layer(numoflayers)%pot(jmi,nn))>maxrk21) then
        !     maxrk21 = real(layer(numoflayers)%pot(jmi,nn))
        ! end if
        ! love_k(jmi) = cmplx(maxrk21 / ff, 0.0_dp, dp)

        jmi = jm(j, 2)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        topur (jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%ur_i(jmi,2)
        ff = real(forcing_k(jm(2,2)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxrk22) then
            maxrk22 = real(layer(numoflayers)%pot(jmi,nn))
            maxrh22 = real(layer(numoflayers)%ur_i(jmi,2))
            maxrl22 = real(layer(numoflayers)%ut_i(jmi,2))
        end if
        ffi = imag(forcing_k(jm(2,2)))
        if(imag(layer(numoflayers)%pot(jmi,nn))>maxik22) then
            maxik22 = imag(layer(numoflayers)%pot(jmi,nn))
            maxih22 = imag(layer(numoflayers)%ur_i(jmi,2))
            maxrl22 = imag(layer(numoflayers)%ut_i(jmi,2))
        end if
        if((ff>0).and.(ffi>0)) then
            love_k(jmi) = cmplx(maxrk22 / ff, maxik22 / ffi, dp)
            love_h(jmi) = cmplx(maxrh22/(ff/layer(numoflayers)%g_i(nn)), maxih22/(ffi/layer(numoflayers)%g_i(nn)), dp)
            love_l(jmi) = cmplx(maxrl22/(ff*sqrt(30.0_dp/5.0_dp)/layer(numoflayers)%g_i(nn)), maxil22/(ff*sqrt(30.0_dp/5.0_dp)/layer(numoflayers)%g_i(nn)), dp)
        end if
    end do
    if(mod(it2, num_time_steps) == 0) then
        maxk20  = -1.0_dp;  maxh20 = -1.0_dp;  maxl20 = -1.0_dp;
        maxrk21 = -1.0_dp; maxrh21 = -1.0_dp; maxrl21 = -1.0_dp;
        maxik21 = -1.0_dp; maxih21 = -1.0_dp; maxil21 = -1.0_dp;
        maxrk22 = -1.0_dp; maxrh22 = -1.0_dp; maxrl22 = -1.0_dp;
        maxik22 = -1.0_dp; maxih22 = -1.0_dp; maxil22 = -1.0_dp;
    end if
end subroutine ind_pot2_coriolis

subroutine ind_pot2(it2)
    implicit none
    integer  :: j, m, it2, nn, jmi
    real(dp) :: ff, ffi

    if(it2==1) then
        maxk20  = -1.0_dp
        maxrk21 = -1.0_dp
        maxik21 = -1.0_dp
        maxrk22 = -1.0_dp
        maxik22 = -1.0_dp
    end if

    j = 2; m = 0
    call induced_potential2(j, m)
    j = 2; m = 1
    call induced_potential2(j, m)
    j = 2; m = 2
    call induced_potential2(j, m)

    nn = layer(numoflayers)%n
    do j = 2, 2
        jmi = jm(j, 0)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff  = real(forcing_k(jm(2,0)))
        ! print*, real(layer(numoflayers)%pot(jmi,nn)), maxk20
        if(real(layer(numoflayers)%pot(jmi,nn))>maxk20) then
            maxk20 = real(layer(numoflayers)%pot(jmi,nn))
            if(ff>0) love_k(jmi) = cmplx(maxk20/ff, 0.0, dp)
        end if

        ! jmi = jm(j, 1)
        ! indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ! ff = real(forcing_k(jm(2,1)))
        ! if(real(layer(numoflayers)%pot(jmi,nn))>maxrk21) then
        !     maxrk21 = real(layer(numoflayers)%pot(jmi,nn))
        ! end if
        ! love_k(jmi) = cmplx(maxrk21 / ff, 0.0_dp, dp)

        jmi = jm(j, 2)
        indpot(jmi, mod(it2-1, num_time_steps) + 1) = layer(numoflayers)%pot(jmi,nn)
        ff = real(forcing_k(jm(2,2)))
        if(real(layer(numoflayers)%pot(jmi,nn))>maxrk22) then
            maxrk22 = real(layer(numoflayers)%pot(jmi,nn))
        end if
        ffi = imag(forcing_k(jm(2,2)))
        if(imag(layer(numoflayers)%pot(jmi,nn))>maxik22) then
            maxik22 = imag(layer(numoflayers)%pot(jmi,nn))
        end if
        if((ff>0).and.(ffi>0))love_k(jmi) = cmplx(maxrk22 / ff, maxik22 / ffi, dp)
    end do
    if(mod(it2, num_time_steps) == 0) then
        maxk20  = -1.0_dp
        maxrk21 = -1.0_dp
        maxik21 = -1.0_dp
        maxrk22 = -1.0_dp
        maxik22 = -1.0_dp
    end if
end subroutine ind_pot2

subroutine induced_potential2(j,m)
    implicit none
    integer     :: il, ir, j, m, jms, irs
    real(dp)    :: r, cg, cj, dr, cint0, cint1, cext0, cext1, scale
    complex(dp) :: above, below, boundary_a, boundary_b
    above = dzero; below = dzero
    jms = jm(j,m)
    cg = 4.0_dp * Grav_Const * Pi
    cj = 1.0_dp / (2.0_dp*real(j,dp)+1.0_dp)

    r  = layer(1)%r_i(1)
    call integrator(r, j, m, 1, 1, above, below)
    call boundary_contribution(r, j, m, 1, boundary_a, boundary_b)
    layer(1)%pot(jms, 1) = (above + boundary_a) * cj * cg * r

    do il = 1, numoflayers
        if(il>1)then
            r  = layer(il)%r_i(1)
            call boundary_contribution(r, j, m, il, boundary_a, boundary_b)
            call integrator(r, j, m, il, 1, above, below)
            layer(il)%pot(jms, 1) = (above + below) * cg * cj * r + (boundary_a + boundary_b) * cj * cg * r
            irs = 2
        else
            irs = 3
        end if
        do ir = irs, layer(il)%n
            r  = layer(il)%r_i(ir)
            dr = layer(il)%r_i(ir) - layer(il)%r_i(ir-1)
            cext1 = (r/layer(il)%r_i(ir  )) ** (j-1)
            cext0 = (r/layer(il)%r_i(ir-1)) ** (j-1)
            scale = (r/layer(il)%r_i(ir-1)) ** (j-1)
            above = above*scale - 0.5_dp*dr*(layer(il)%rhoe(ir, jms)*cext1 + layer(il)%rhoe(ir-1, jms)*cext0)
            cint1 = (layer(il)%r_i(ir)  /r)  ** (j+2)
            cint0 = (layer(il)%r_i(ir-1)/r)  ** (j+2)
            scale = (layer(il)%r_i(ir-1)/r)  ** (j+2)
            below = below*scale + 0.5_dp*dr*(layer(il)%rhoe(ir, jms)*cint1 + layer(il)%rhoe(ir-1, jms)*cint0)
            call boundary_contribution(r, j, m, il, boundary_a, boundary_b)
            layer(il)%pot(jms, ir) = (above + below) * cg * cj * r + (boundary_a + boundary_b) * cj * cg * r
        end do
    end do
end subroutine induced_potential2

subroutine boundary_contribution(r, j, m, layer_number, ind_pot_a, ind_pot_b)
    implicit none
    integer     :: j, m, il, layer_number, nn, jms
    real(dp)    :: drho, r
    complex(dp) :: ind_pot_b, ind_pot_a
    ind_pot_a = dzero
    ind_pot_b = dzero
    jms = jm(j,m)
    if(layer_number > 1) then
        do il = 1, layer_number-1
            nn      = layer(il)%n
            drho    = layer(il)%rho_i(nn) - layer(il+1)%rho_i(1)
            ind_pot_b = ind_pot_b + drho * (layer(il)%ur_i(jms, 2)) * (layer(il)%r_i(nn)  /r)  ** (j+2)
        end do
    end if
    do il = layer_number, numoflayers
        nn = layer(il)%n
        if((il<numoflayers)) then
            drho    = layer(il)%rho_i(nn) - layer(il+1)%rho_i(1)
        elseif(il==numoflayers)then
            drho    = layer(il)%rho_i(nn)
        end if
        ind_pot_a = ind_pot_a + drho * (layer(il)%ur_i(jms, 2)) * (r/layer(il)%r_i(nn)  )** (j-1)
    end do
end subroutine boundary_contribution

subroutine integrator(r, j, m, layer_number, radial_number, ind_pot_a, ind_pot_b)
    implicit none
    integer    :: layer_number, radial_number, il, ir, j, m, jms, nr
    real(dp)    :: cint0, cint1, cext0, cext1, r, dr
    complex(dp) :: ind_pot_a, ind_pot_b

    !! Avoid the integral in incompressible models
    if (SWITCH_compressible == 0) integral_on = 0
    if (SWITCH_compressible == 1) integral_on = 1
    integral_on = 1
    !! INTEGRATION OF THE PARTS THAT ARE BELOW THE GIVEN POINT
    jms = jm(j, m)
    ind_pot_b = 0.0_dp
    if(integral_on == 1) then
        do il = 1, layer_number
            if(layer_number > 1) then
                if(il < layer_number) then
                    nr = layer(il)%n
                else 
                    nr = radial_number
                end if
            else
                nr = radial_number
            end if
            do ir = 1, nr-1
                dr = layer(il)%r_i(ir+1) - layer(il)%r_i(ir)
                cint0 = (layer(il)%r_i(ir)  /r)  ** (j+2)
                cint1 = (layer(il)%r_i(ir+1)/r)  ** (j+2)
                ind_pot_b = ind_pot_b + 0.5_dp * dr*  (cint0 * layer(il)%rhoe(ir, jms) + cint1 * layer(il)%rhoe(ir+1, jms))
            end do
        end do
    end if
    !!Boundary term
    ! if(layer_number > 1) then
    !     do il = 1, layer_number-1
    !         nn      = layer(il)%n
    !         drho    = layer(il)%rho_i(nn) - layer(il+1)%rho_i(1)
    !         ind_pot = ind_pot + drho * (layer(il)%ur_i(jms, 2)) * (layer(il)%r_i(nn)  /r)  ** (j+2)
    !     end do
    ! end if

    !! INTEGRATION OF THE PARTS THAT ARE ABOVE THE GIVEN POINT
    ind_pot_a = dzero
    if(integral_on == 1) then
        if(layer_number>=2)then
            do il = layer_number, numoflayers
                if(layer_number < numoflayers) then
                    if(il > layer_number) then
                        nr = 1
                    else
                        nr = radial_number
                    end if
                else
                    nr = radial_number
                end if
                do ir = nr, layer(il)%n-1
                    dr = layer(il)%r_i(ir+1) - layer(il)%r_i(ir)
                    cext0 = (r/layer(il)%r_i(ir)  )** (j-1)
                    cext1 = (r/layer(il)%r_i(ir+1))** (j-1)
                    ind_pot_a = ind_pot_a + 0.5_dp * dr * (cext0 * layer(il)%rhoe(ir, jms) + cext1 * layer(il)%rhoe(ir+1, jms))
                end do
            end do
        end if
    end if
    ! do il = layer_number, numoflayers
    !     nn = layer(il)%n
    !     if((il<numoflayers)) then
    !         drho    = layer(il)%rho_i(nn) - layer(il+1)%rho_i(1)
    !     elseif(il==numoflayers)then
    !         drho    = layer(il)%rho_i(nn)
    !     end if
    !     ind_pot = ind_pot + drho * (layer(il)%ur_i(jms, 2)) * (r/layer(il)%r_i(nn)  )** (j-1)
    ! end do
end subroutine integrator


subroutine love_numbers
    implicit none
    integer    :: nn, j
    complex(dp) :: surface_pot

    nn = layer(numoflayers)%n

    surface_pot = layer(numoflayers)%pot(nn, jm2(2,0))
    k20 = real(surface_pot)

    surface_pot = layer(numoflayers)%pot(nn, jm2(2,2))
    real_k22    = real(surface_pot)
    imag_k22    = imag(surface_pot)

    h20         = real(layer(numoflayers)%ur_i(nn, jm2(2,0)))!/real(forcing_pot/layer(numoflayers)%g_i(nn))
    real_h22    = real(layer(numoflayers)%ur_i(nn, jm2(2,2)))
    imag_h22    = imag(layer(numoflayers)%ur_i(nn, jm2(2,2)))

    if(SWITCH_velocity==1)then
    l20         = real(layer(numoflayers)%dis(nn, jml2(2,0,1))*sqrt(3.0_dp/5.0_dp) + layer(numoflayers)%dis(nn, jml2(2,0,3))*sqrt(2.0_dp/5.0_dp))
    real_l22    = real(layer(numoflayers)%dis(nn, jml2(2,2,1))*sqrt(3.0_dp/5.0_dp) + layer(numoflayers)%dis(nn, jml2(2,2,3))*sqrt(2.0_dp/5.0_dp))
    imag_l22    = imag(layer(numoflayers)%dis(nn, jml2(2,2,1))*sqrt(3.0_dp/5.0_dp) + layer(numoflayers)%dis(nn, jml2(2,2,3))*sqrt(2.0_dp/5.0_dp))
    else
    l20         = real(layer(numoflayers)%vel(nn, jml2(2,0,1))*sqrt(3.0_dp/5.0_dp) + layer(numoflayers)%vel(nn, jml2(2,0,3))*sqrt(2.0_dp/5.0_dp))
    real_l22    = real(layer(numoflayers)%vel(nn, jml2(2,2,1))*sqrt(3.0_dp/5.0_dp) + layer(numoflayers)%vel(nn, jml2(2,2,3))*sqrt(2.0_dp/5.0_dp))
    imag_l22    = imag(layer(numoflayers)%vel(nn, jml2(2,2,1))*sqrt(3.0_dp/5.0_dp) + layer(numoflayers)%vel(nn, jml2(2,2,3))*sqrt(2.0_dp/5.0_dp))
    end if


    do j = 1, 3
        kh(3*(j-1)+1) = real(layer(numoflayers)%pot(nn, jm2(2*(j+1), 0))) !!4,0 
        kh(3*(j-1)+2) = real(layer(numoflayers)%pot(nn, jm2(2*(j+1), 2))) !!4,2, real 
        kh(3*(j-1)+3) = imag(layer(numoflayers)%pot(nn, jm2(2*(j+1), 2))) !!4,2, imag 
    end do
end subroutine love_numbers


complex(dp) function tides(den, rr, m, tt, gg) result (tid)
    implicit none
    integer          :: m, it2
    real(dp)         :: c20, rec22, imc22, den, rr, tt, rec21, w2
    character(len=1) :: gg   !!g for gradient and p for potential
    w2 = ang_vel * ang_vel
    it2 = nint(tt / time_step)
    if(ecc > 0.0_dp)then
        ecc = ecc
    else 
        ecc = ecc_arr(it2)
    end if
    ! print*, tt - time_ecc(it2), ecc
    if(gg == 'g')then
        c20    =  w2 * ecc * sqrt(18.0_dp * Pi) * cos(ang_vel * tt)
        rec21  =  w2 * obl * sqrt(12.0_dp * Pi) * cos(ang_vel * tt)
        rec22  =  w2 * ecc * sqrt(27.0_dp * Pi) * cos(ang_vel * tt)
        imc22  =  w2 * ecc * sqrt(48.0_dp * Pi) * sin(ang_vel * tt)
    else
        c20    =  (rr/den) * w2 * ecc * sqrt(18.0_dp * Pi / 10.0_dp) * cos(ang_vel * tt)
        rec21  =  (rr/den) * w2 * obl * sqrt(12.0_dp * Pi / 10.0_dp) * cos(ang_vel * tt)
        rec22  =  (rr/den) * w2 * ecc * sqrt(27.0_dp * Pi / 10.0_dp) * cos(ang_vel * tt)
        imc22  =  (rr/den) * w2 * ecc * sqrt(48.0_dp * Pi / 10.0_dp) * sin(ang_vel * tt)
    end if
    
    tid = dzero
    if    (m.eq.0)then
        tid =  den*rr*cmplx(c20,0.0, dp)
    elseif(m.eq.1)then
        tid =  den*rr*cmplx(rec21, 0.0_dp,dp)
    elseif(m.eq.2)then
        tid = -den*rr*cmplx(rec22,-imc22, dp)
    end if
end function tides

complex(dp) function mm_tides(den, rr, m, tt, gg) result (tid)
    implicit none
    integer    :: m
    real(dp)    :: den, rr, tt
    character(len=1) :: gg !!g for gradient and p for potential
    
    Europa%mass    = 4.79984e22_dp
    Europa%period  = rotation_period*0.5_dp
    Europa%ang_vel = 2.0_dp*pi/Europa%period
    Europa%semimajor_a = 671100.0e3_dp
    Io%mass        = 8.93e22_dp
    Io%period      = rotation_period*0.25_dp
    Io%ang_vel     = 2.0_dp*pi/Io%period
    Io%semimajor_a = 421700.e3_dp
    Europa%nij = -(Europa%ang_vel-ang_vel)
    Io%nij     = -(Io%ang_vel-ang_vel)
    Europa%rij = sqrt(Europa%semimajor_a**2 + semimajor_a**2 - 2.0_dp*Europa%semimajor_a*semimajor_a*cos(Europa%nij*tt))
    Io%rij     = sqrt(Io%semimajor_a**2     + semimajor_a**2 - 2.0_dp*Io%semimajor_a    *semimajor_a*cos(Io%nij*tt - Pi))!!with initial delay of Pi

    Europa%cospij = (semimajor_a - Europa%semimajor_a * cos(Europa%nij * tt)) / Europa%rij
    Europa%sinpij =                Europa%semimajor_a * sin(Europa%nij * tt)  / Europa%rij

    Io%cospij     = (semimajor_a - Io%semimajor_a * cos(Io%nij * tt - Pi)) / Io%rij !!with initial delay of Pi
    Io%sinpij     =                Io%semimajor_a * sin(Io%nij * tt - Pi)  / Io%rij !!with initial delay of Pi

    Europa%cos2pij = Europa%cospij**2 - Europa%sinpij**2
    Europa%sin2pij = Europa%sinpij*Europa%cospij*2.0_dp

    Io%cos2pij = Io%cospij**2 - Io%sinpij**2
    Io%sin2pij = Io%sinpij*Io%cospij*2.0_dp

    !!Compute the gradient...
    if(gg == 'g')then
        Europa%C20 = sqrt(2.0_dp * Pi) * Grav_Const * Europa%mass / (Europa%rij**3) !sqrt(pi/5.0_dp)*Grav_Const*Europa%mass/(Europa%rij**3)
        Europa%C22 = sqrt(3.0_dp * Pi) * Grav_Const * Europa%mass * Europa%cos2pij / (Europa%rij**3)
        Europa%S22 = sqrt(3.0_dp * Pi) * Grav_Const * Europa%mass * Europa%sin2pij / (Europa%rij**3)

        Io%C20     = sqrt(2.0_dp * Pi) * Grav_Const * Io%mass     / (Io%rij**3) 
        Io%C22     = sqrt(3.0_dp * Pi) * Grav_Const * Io%mass * Io%cos2pij / (Io%rij**3)
        Io%S22     = sqrt(3.0_dp * Pi) * Grav_Const * Io%mass * Io%sin2pij / (Io%rij**3)
    else
        Europa%C20 = (rr/den) * sqrt(        Pi /  5.0_dp) * Grav_Const * Europa%mass                  / (Europa%rij**3)
        Europa%C22 = (rr/den) * sqrt(3.0_dp * Pi / 10.0_dp) * Grav_Const * Europa%mass * Europa%cos2pij / (Europa%rij**3)
        Europa%S22 = (rr/den) * sqrt(3.0_dp * Pi / 10.0_dp) * Grav_Const * Europa%mass * Europa%sin2pij / (Europa%rij**3)

        Io%C20     = (rr/den) * sqrt(        Pi /  5.0_dp) * Grav_Const * Io%mass                      / (Io%rij**3)
        Io%C22     = (rr/den) * sqrt(3.0_dp * Pi / 10.0_dp) * Grav_Const * Io%mass * Io%cos2pij         / (Io%rij**3)
        Io%S22     = (rr/den) * sqrt(3.0_dp * Pi / 10.0_dp) * Grav_Const * Io%mass * Io%sin2pij         / (Io%rij**3)
    end if

    tid = dzero
    if    (m.eq.0)then
        tid =  ( den*rr*cmplx(Europa%C20+Io%C20, 0.0_dp, dp))
    elseif(m.eq.1)then
    else
        tid =  (-den*rr*cmplx(Europa%C22+Io%C22, -Europa%S22-Io%S22, dp))
    end if
end function mm_tides


complex(dp) function lunar_tides(den, rr, m, tt, gg) result(lun_tid)
    implicit none
    integer :: m
    real(dp) :: den, tt, rr
    real(dp) :: C22, S22
    character(len=1) :: gg !!g for gradient and p for potential
    
    lun_tid = dzero
    Moon%mass   = 7.3e22_dp
    Moon%rij    = 0.38e5_dp * km_to_m !!Distance to the Moon in meters -- minimum 0.38 and maximum 38.0 -- range taken from Katz et al. 2025
    Moon%nij    = 7.0e-5_dp           !! Earth rotation -- minimum 7.0 and maximum 70 -- range taken from Katz et al. 2025
    ! C22 = 0.25_dp * Grav_Const * Moon%mass * sqrt(24.0_dp * Pi / 5.0_dp) * cos(2.0_dp * Moon%nij * tt)/ (Moon%rij**3)
    C22 = 0.25_dp * Grav_Const * Moon%mass * sqrt(24.0_dp * Pi / 5.0_dp) / (Moon%rij**3)
    S22 = 0.25_dp * Grav_Const * Moon%mass * sqrt(24.0_dp * Pi / 5.0_dp) * sin(2.0_dp * Moon%nij * tt)/ (Moon%rij**3)

    if(gg == 'g') then
        if(m==0)then
        else if (m==1)then
        else
            lun_tid = den*rr*C22* 5.0_dp * sqrt(0.4_dp) * exp(-2.0_dp * im * Moon%nij * tt)
        end if
    else
        if(m==0)then
        else if (m==1)then
        else
            lun_tid = (rr**2)*C22*exp(-2.0_dp * im * Moon%nij * tt)
        end if
    end if
end function lunar_tides

end module mod_potential