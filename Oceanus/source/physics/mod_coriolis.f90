module mod_coriolis
    use mod_parameters
    use mod_variables
    use mod_math
    use mod_SHTns
    implicit none
    
contains

    subroutine convert2jmax
        implicit none
        
        if(SWITCH_coriolis == 1) then
            jstart=1; jfinal=jmax
        end if
        status = 1
    end subroutine convert2jmax


    subroutine coriolis3_f
        implicit none
        integer  :: jmind, il, compute_vector, oceindx
        real(dp) :: j, m, ff
        complex(dp) :: fcjm1_vjm1, fcjm1_vj
        complex(dp) :: fcj_vjm1, fcj_vj, fcj_vjp1
        complex(dp) :: fcjp1_vj, fcjp1_vjp1
        do il = 1, numoflayers
            layer(il)%fc(:,:) = dzero
        end do

        ff = -2.0_dp *  ang_vel
        compute_vector = jml(jmax, jmax, jmax+1)
        !$OMP PARALLEL
        do il = 1, num_oceans
            oceindx = layer_ocean(il)
            !$OMP DO PRIVATE(j, m, fcjm1_vjm1, fcjm1_vj, fcj_vjm1, fcj_vj, fcj_vjp1, fcjp1_vj, fcjp1_vjp1)
            do jmind = 2, compute_vector, 3
                j = real(jmlindx(jmind,2), dp); m = real(jmlindx(jmind,3), dp)
                fcjm1_vjm1 =     (Im/j) * sqrt( (j-1.0_dp) * (j**2 - m**2) / (2.0_dp*j-1.0_dp) )
                fcjm1_vj   = -(Im*m/j) 
                
                fcj_vjm1   =     (Im/j)     * sqrt( (j+1.0_dp) * (j**2 - m**2) / (2.0_dp*j+1.0_dp) )
                fcj_vj     = -(Im*m/(j*(j+1.0_dp)))
                fcj_vjp1   = (Im/(j+1.0_dp)) * sqrt( (j*((j+1.0_dp)**2 - m**2)) / (2.0_dp*j+1.0_dp) )
                
                fcjp1_vj   = (Im*m)/(j+1.0_dp)
                fcjp1_vjp1 = (Im/(j+1.0_dp)) * sqrt( ((j+2.0_dp)*((j+1.0_dp)**2 - m**2)) / (2.0_dp*j+3.0_dp) )

                layer(il)%fc(jmind, :) = layer(il)%rho_l(:) * ff * (fcjm1_vj   * layer(il)%vel(jmind, :))
                if(int(j)> 1) layer(il)%fc(jmind, :) = layer(il)%fc(jmind, :) + layer(il)%rho_l(:) * ff * (fcjm1_vjm1 * layer(il)%vel(jml(int(j)-1, int(m), int(j)-1), :))

                layer(il)%fc(jmind+1, :) = layer(il)%rho_l(:) * ff * (fcj_vj     * layer(il)%vel(jmind+1, :) )
                if(int(j)<jmax) layer(il)%fc(jmind+1, :) = layer(il)%fc(jmind+1, :) + layer(il)%rho_l(:) * ff * fcj_vjp1 * layer(il)%vel(jml(int(j)+1, int(m), int(j)), :)
                if(int(j)>1   ) layer(il)%fc(jmind+1, :) = layer(il)%fc(jmind+1, :) + layer(il)%rho_l(:) * ff * fcj_vjm1 * layer(il)%vel(jml(int(j)-1, int(m), int(j)), :)
                
                layer(il)%fc(jmind+2, :) = layer(il)%rho_l(:) * ff * (fcjp1_vj   * layer(il)%vel(jmind+2, :) ) 
                if(int(j)<jmax) layer(il)%fc(jmind+2, :) = layer(il)%fc(jmind+2, :) + layer(il)%rho_l(:) * ff * fcjp1_vjp1 * layer(il)%vel(jml(int(j)+1 ,int(m),int(j)+1), :)
            end do
            !$OMP END DO
        end do
        !$OMP END PARALLEL
    end subroutine coriolis3_f

    subroutine coriolis3_r
        implicit none
        integer  :: jmind, il, oceindx
        real(dp) :: j, m, ff
        complex(dp) :: fcjm1_vjm1, fcjm1_vj
        complex(dp) :: fcj_vjm1, fcj_vj, fcj_vjp1
        complex(dp) :: fcjp1_vj, fcjp1_vjp1
        do il = 1, numoflayers
            layer(il)%fc(:,:) = dzero
        end do
        ff = -2.0_dp *  ang_vel
        
        !$OMP PARALLEL
        do il = 1, num_oceans
            oceindx = layer_ocean(il)
            
            !$OMP DO PRIVATE(j, m, fcjm1_vjm1, fcjm1_vj, fcj_vjm1, fcj_vj, fcj_vjp1, fcjp1_vj, fcjp1_vjp1)
            do jmind = 2, compute_vector_harm, 3
                j = real(jmlindx(jmind,2), dp); m = real(jmlindx(jmind,3), dp)
                fcjm1_vjm1 =     (Im/j) * sqrt( (j-1.0_dp) * (j**2 - m**2) / (2.0_dp*j-1.0_dp) )
                fcjm1_vj   = -(Im*m/j) 
                
                fcj_vjm1   =     (Im/j)     * sqrt( (j+1.0_dp) * (j**2 - m**2) / (2.0_dp*j+1.0_dp) )
                fcj_vj     = -(Im*m/(j*(j+1.0_dp)))
                fcj_vjp1   = (Im/(j+1.0_dp)) * sqrt( (j*((j+1.0_dp)**2 - m**2)) / (2.0_dp*j+1.0_dp) )
                
                fcjp1_vj   = (Im*m)/(j+1.0_dp)
                fcjp1_vjp1 = (Im/(j+1.0_dp)) * sqrt( ((j+2.0_dp)*((j+1.0_dp)**2 - m**2)) / (2.0_dp*j+3.0_dp) )
                ! print*, j, m, fcj_vjm1, fcj_vj, fcj_vjp1
                layer(oceindx)%fc(jmind, :) = layer(oceindx)%rho_l(:) * ff * (fcjm1_vj   * layer(oceindx)%vel(jmind, :))
                if(int(j)> 1) layer(oceindx)%fc(jmind, :) = layer(oceindx)%fc(jmind, :) + layer(oceindx)%rho_l(:) * ff * (fcjm1_vjm1 * layer(oceindx)%vel(jml(int(j)-1, int(m), int(j)-1), :))

                layer(oceindx)%fc(jmind+1, :) = layer(oceindx)%rho_l(:) * ff * (fcj_vj     * layer(oceindx)%vel(jmind+1, :) )
                if(int(j)<jmax) layer(oceindx)%fc(jmind+1, :) = layer(oceindx)%fc(jmind+1, :) + layer(oceindx)%rho_l(:) * ff * fcj_vjp1 * layer(oceindx)%vel(jml(int(j)+1, int(m), int(j)), :)
                if(int(j)>1   ) layer(oceindx)%fc(jmind+1, :) = layer(oceindx)%fc(jmind+1, :) + layer(oceindx)%rho_l(:) * ff * fcj_vjm1 * layer(oceindx)%vel(jml(int(j)-1, int(m), int(j)), :)
                
                layer(oceindx)%fc(jmind+2, :) = layer(oceindx)%rho_l(:) * ff * (fcjp1_vj   * layer(oceindx)%vel(jmind+2, :) ) 
                if(int(j)<jmax) layer(oceindx)%fc(jmind+2, :) = layer(oceindx)%fc(jmind+2, :) + layer(oceindx)%rho_l(:) * ff * fcjp1_vjp1 * layer(oceindx)%vel(jml(int(j)+1 ,int(m),int(j)+1), :)
            end do
            !$OMP END DO
        end do
        !$OMP END PARALLEL
    end subroutine coriolis3_r

    subroutine coriolis2
        implicit none
        integer ::  j,  m, il, jm1, jm0, jp1
        real(dp) :: jj, mm, ff
        complex(dp) :: fcjm1_vjm1, fcjm1_vj
        complex(dp) :: fcj_vjm1, fcj_vj, fcj_vjp1
        complex(dp) :: fcjp1_vj, fcjp1_vjp1
        do il = 1, numoflayers
            layer(il)%fc(:,:) = dzero
        end do
        
        ff = -2.0_dp *  ang_vel
        do j = 1, jmax
            jj = real(j, dp)
            do m = 0, j
                mm = real(m, dp)
                jm1=jml(j,m, j-1); jm0=jml(j,m,j); jp1=jml(j,m,j+1)
                fcjm1_vjm1 =     (Im/jj) * sqrt( (jj-1.0_dp) * (jj**2 - mm**2) / (2.0_dp*jj-1.0_dp) )
                fcjm1_vj   = -(Im*mm/jj) 

                fcj_vjm1   =     (Im/jj)     * sqrt( (jj+1.0_dp) * (jj**2 - mm**2) / (2.0_dp*jj+1.0_dp) )
                fcj_vj     = -(Im*mm/(jj*(jj+1.0_dp)))
                fcj_vjp1   = (Im/(jj+1.0_dp)) * sqrt( (jj*((jj+1.0_dp)**2 - mm**2)) / (2.0_dp*jj+1.0_dp) )

                fcjp1_vj   = (Im*mm)/(jj+1.0_dp)
                fcjp1_vjp1 = (Im/(jj+1.0_dp)) * sqrt( ((jj+2.0_dp)*((jj+1.0_dp)**2 - mm**2)) / (2.0_dp*jj+3.0_dp) )

                do il = 1, 2
                    layer(il)%fc(:, jm1) = layer(il)%rho_l(:) * ff * (fcjm1_vj   * layer(il)%vel(:, jml(j   ,m,j-1)))
                    if(j> 1) layer(il)%fc(:, jm1) = layer(il)%fc(:, jm1) + layer(il)%rho_l(:) * ff * (fcjm1_vjm1 * layer(il)%vel(:, jml(j-1,m, j-1)))

                    layer(il)%fc(:, jm0) = layer(il)%rho_l(:) * ff * (fcj_vj     * layer(il)%vel(:, jml(j   ,m,j  )) )
                    if(j<jmax) layer(il)%fc(:, jm0) = layer(il)%fc(:, jm0) + layer(il)%rho_l(:) * ff * fcj_vjp1 * layer(il)%vel(:, jml(j+1, m, j))
                    if(j>1   ) layer(il)%fc(:, jm0) = layer(il)%fc(:, jm0) + layer(il)%rho_l(:) * ff * fcj_vjm1 * layer(il)%vel(:, jml(j-1, m, j))
                    
                    layer(il)%fc(:, jp1) = layer(il)%rho_l(:) * ff * (fcjp1_vj   * layer(il)%vel(:, jml(j  ,m, j+1)) ) 
                    if(j<jmax) layer(il)%fc(:, jp1) = layer(il)%fc(:, jp1) + layer(il)%rho_l(:) * ff * fcjp1_vjp1 * layer(il)%vel(:, jml(j+1 ,m,j+1))

                end do
            end do
        end do

    end subroutine coriolis2

    subroutine coriolis
        implicit none
        real(dp) :: cc, c6, co, coef1, fco
        integer :: j2, deg, ord, l, jmli, jml1, i, k
        do i = 1, numoflayers
            layer(i)%fc(:,:) = dcmplx(0.0_dp, 0.0_dp)
        end do

        do deg =  1, jmax
            order: do ord = 0, min(2, deg), 2
                vec_l: do l = abs(deg - 1), deg + 1
                    jmli   = jml(deg, ord, l)
                    coef1 = (-1.0_dp)**(deg + l)
                    ddeg_2: do j2 = max(deg-1, 1, ord, l-1), min(deg+1, jmax, l+1)
                    jml1 = jml(j2, ord, l)
                    call cleb_1(j2, ord, 1, 0, deg, ord, cc)
                    call six_j(j2, 1, l, deg, 1, c6)
                    do i = 1, 2
                        do k = 1, layer(i)%n+1
                            co = 2.0_dp * sqrt(6.0_dp) * ang_vel * layer(i)%rho_l(k)
                            fco = co * coef1 * sqrt(real(2*j2 + 1, dp)) * cc * c6
                            layer(i)%fc(jmli, k) = layer(i)%fc(jmli, k) + (Im * fco * layer(i)%vel(jml1, k))
                        end do
                    end do
                end do ddeg_2
                end do vec_l
            end do order
        end do
    end subroutine coriolis


    subroutine coriolis_pvel
        implicit none
        real(dp) :: cc, c6, co, coef1, fco
        integer :: j2, deg, ord, l, jmlc, jml1, i, k
        do i = 1, numoflayers
            layer(i)%fc0(:,:) = dcmplx(0.0_dp, 0.0_dp)
        end do

        do deg =  0, jmax
            order: do ord = 0, min(2,deg), 2
                vec_l: do l = abs(deg - 1), deg + 1
                    jmlc   = jml2(deg, ord, l)
                    coef1 = (-1.0_dp)**(deg + l)
                    ddeg_2: do j2 = max(deg-1, 1, ord, l-1), min(deg+1, jmax, l+1)
                        jml1 = jml2(j2, ord, l)
                        call cleb_1(j2, ord, 1, 0, deg, ord, cc)
                        call six_j(j2, 1, l, deg, 1, c6)
                        do i = 1, 2
                            do k = 1, layer(i)%n-1
                                co = 2.0_dp * sqrt(6.0_dp) * ang_vel * layer(i)%rho_l(k+1)
                                fco = co * coef1 * sqrt(real(2*j2 + 1, dp)) * cc * c6
                                layer(i)%fc0(k, jmlc) = layer(i)%fc0(k, jmlc) + (Im * fco * layer(i)%pvel(k+1, jml1))
                            end do
                        end do
                    end do ddeg_2
                end do vec_l
            end do order
        end do
    end subroutine coriolis_pvel



    subroutine six_j(j1, j2, j, l2, l, c6)
        !     6-j symbol  :  {j1,j2,j
        !                     1,l2,l}
        !.....................................................................
        implicit none
        integer :: j1, j2, j, l2, l
        real(dp) :: c6, a, b, c, s, zn
        real(dp) :: hh, dd

        c6=0.0_dp
    
        if(l2.gt.(j+1).or.iabs(j-1).gt.l2)   return
        if(j.gt.(j1+j2).or.iabs(j1-j2).gt.j) return
        if(l.gt.(j2+1).or.iabs(j2-1).gt.l)   return
    
        a=real(j1,dp)
        b=real(j2,dp)
        c=real(j,dp)
        s=a+b+c
        zn=(-1.0_dp)**(j1+j2+j)
    
        if((j-l2).lt.0) then
            if((j2-l).lt.0)then
                hh=(s+2)*(s+3)*(s-a-a+1)*(s-a-a+2)
                dd=(b+b+1)*(b+1)*(b+b+3)*(c+c+1)*(c+1)*(c+c+3)
            end if
            if((j2-l).eq.0)then
                hh=(s+2)*(s-c-c)*(s-b-b+1)*(s-a-a+1)
                dd=b*(b+b+1)*(b+1)*(c+c+1)*(c+1)*(c+c+3)
                zn=-zn
            end if
            if((j2-l).gt.0)then
                hh=(s-c-c-1)*(s-c-c)*(s-b-b+1)*(s-b-b+2)
                dd=(b+b-1)*b*(b+b+1)*(c+c+1)*(c+1)*(c+c+3)
            end if
        end if
    
        if((j-l2).eq.0) then
            if((j2-l).lt.0)then
                hh=(s+2)*(s-c-c+1)*(s-b-b)*(s-a-a+1)
                dd=(b+b+1)*(b+1)*(b+b+3)*c*(c+c+1)*(c+1)
                zn=-zn
            end if
            if((j2-l).eq.0)then
                hh=-a*(a+1)+b*(b+1)+c*(c+1)
                dd=b*(b+b+1)*(b+1)*c*(c+c+1)*(c+1)
                c6=-zn*hh/sqrt(real(dd, dp))/2._dp; return
            end if
            if((j2-l).gt.0)then
                hh=(s+1)*(s-c-c)*(s-b-b+1)*(s-a-a)
                dd=(b+b-1)*b*(b+b+1)*c*(c+c+1)*(c+1)
            end if
        end if
    
        if((j-l2).gt.0) then
            if((j2-l).lt.0)then
                hh=(s-c-c+1)*(s-c-c+2)*(s-b-b-1)*(s-b-b)
                dd=(b+b+1)*(b+1)*(b+b+3)*(c+c-1)*c*(c+c+1)
            end if
            if((j2-l).eq.0)then
                hh=(s+1)*(s-c-c+1)*(s-b-b)*(s-a-a)
                dd=b*(b+b+1)*(b+1)*(c+c-1)*c*(c+c+1)
            end if
            if((j2-l).gt.0)then
                hh=s*(s+1)*(s-a-a-1)*(s-a-a)
                dd=(b+b-1)*b*(b+b+1)*(c+c-1)*c*(c+c+1)
            end if
        end if
    
        c6=zn*sqrt((real(hh, dp))/real(dd, dp))*0.5_dp
    end subroutine six_j

    subroutine cleb_1(j1,m1,j2,m2,j,m,cg)
        !
        !      j m
        !     C
        !      j1 m1 1 m2
        !........................................................................
        implicit none
        integer        :: j1, m1, j2, m2, j, m
        real(dp) :: cg, c, g, zn

        cg=0.0_dp
        if(j2.ne.1) return
        if(abs(j1-j).gt.1.or.(j1+j).eq.0) return
        if(abs(m2).gt.1.or.abs(m1).gt.j1) return

        if(m.ne.(m1+m2)) return
        c=real(j, dp)
        g=real(m, dp)
        zn=1.0_dp
        
        if(m2.lt.0)then
            if((j1-j).lt.0)then
                cg=(c-g-1.0_dp)*(c-g)/(c+c-1.0_dp)/(c+c)
            elseif((j1-j).eq.0)then
                cg=(c+g+1.0_dp)*(c-g)/(c+1.0_dp)/(c+c)
            else
                cg=(c+g+2.0_dp)*(c+g+1.0_dp)/(c+c+2.0_dp)/(c+c+3.0_dp)
            end if
        end if

        if(m2.eq.0)then
            if((j1-j).lt.0)then
                cg=(c+g)*(c-g)/(c+c-1.0_dp)/c
            elseif((j1-j).eq.0)then
                cg=g/sqrt(c*(c+1.0_dp)); return
            else
                cg=(c+g+1.0_dp)*(c-g+1.0_dp)/(c+1.0_dp)/(c+c+3.0_dp)
                zn=-1.0_dp
            end if
        end if

        if(m2.gt.0)then
            if((j1-j).lt.0)then
                cg=(c+g-1.0_dp)*(c+g)/(c+c-1.0_dp)/(c+c)
            elseif((j1-j).eq.0)then
                cg=(c+g)*(c-g+1.0_dp)/(c+1.0_dp)/(c+c)
                zn=-1.0_dp
            else
                cg=(c-g+1.0_dp)*(c-g+2.0_dp)/(c+c+2.0_dp)/(c+c+3.0_dp)
            end if
        end if

        cg=zn*sqrt(cg)
    end subroutine cleb_1
    
end module mod_coriolis