module mod_LU
    use mod_parameters
    implicit none
    contains
    
    SUBROUTINE lu_decomp(a,n,m1,m2,np,mp,al,mpl,indx,d)
        implicit none
        INTEGER :: m1,m2,mp,mpl,n,np,indx(n)
        real(dp) :: d,a(np,mp),al(np,mpl)
        real(dp), PARAMETER :: TINY=1.e-30_dp
        INTEGER ::  i,j,k,l,mm
        real(dp) :: dum
        mm=m1+m2+1
        if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) print *,'bad args in bandec'
        l=m1
        do i=1,m1
            do j=m1+2-i,mm
                a(i,j-l)=a(i,j)
            enddo
            l=l-1
            do j=mm-l,mm
                a(i,j)=0.d0
            enddo
        enddo
        d=1.d0
        l=m1
        do k=1,n
            dum=a(k,1)
            i=k
            if(l.lt.n)l=l+1
            do j=k+1,l
                if(abs(a(j,1)).gt.abs(dum))then
                    dum=a(j,1)
                    i=j
                endif
            enddo
            indx(k)=i
            ! if(dum.eq.0.0_dp) a(k,1)=TINY
            if (abs(dum)<1e-16_dp) a(k,1)=TINY
            ! if(dum <= TINY) a(k,1)=TINY
            if(i.ne.k)then
                d=-d
                do j=1,mm
                    dum=a(k,j)
                    a(k,j)=a(i,j)
                    a(i,j)=dum
                enddo
            endif
            do i=k+1,l
                dum=a(i,1)/a(k,1)
                al(k,i-k)=dum
                do j=2,mm
                    a(i,j-1)=a(i,j)-dum*a(k,j)
                enddo
                a(i,mm)=0.d0
            enddo
        enddo
    END SUBROUTINE
    
    SUBROUTINE lu_solve(a,n,m1,m2,np,mp,al,mpl,indx,b)
        implicit none
        INTEGER  :: m1,m2,mp,mpl,n,np,indx(n)
        real(dp) :: a(np,mp),al(np,mpl)
        complex(dp) :: b(n)
        INTEGER :: i,k,l,mm
        complex(dp) :: dum
        mm=m1+m2+1
        if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) print *,'bad args in banbks'
        l=m1
        do k=1,n
            i=indx(k)
            if(i.ne.k)then
                dum=b(k)
                b(k)=b(i)
                b(i)=dum
            endif
            if(l.lt.n)l=l+1
            do i=k+1,l
                b(i)=b(i)-al(k,i-k)*b(k)
            enddo
        enddo
        l=1
        do i=n,1,-1
            dum=b(i)
            do k=2,l
                dum=dum-a(i,k)*b(k+i-1)
            enddo
            b(i)=dum/a(i,1)
            if(l.lt.mm) l=l+1
        enddo
    END SUBROUTINE
    
    SUBROUTINE banmul(a,n,m1,m2,np,mp,x,b)
        implicit none
        INTEGER m1,m2,mp,n,np
        real(dp) a(np,mp),b(n),x(n)
        INTEGER i,j,k
        do i=1,n
        b(i)=0.d0
        k=i-m1-1
        do j=max(1,1-k),min(m1+m2+1,n-k)
            b(i)=b(i)+a(i,j)*x(j+k)
        enddo
        enddo
    END SUBROUTINE
    
    ! SUBROUTINE evalRcond(a,n,m1,m2,rcond)
    ! real(8) :: a(:,:),rcond
    ! integer n,m1,m2
    ! intent(in) :: a,n,m1,m2
    ! intent(out) :: rcond
    ! rcond=0.
    ! END SUBROUTINE
    
    end module
    