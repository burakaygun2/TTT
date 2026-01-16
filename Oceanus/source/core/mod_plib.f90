module mod_plib
    use mod_parameters
    implicit none
    
contains
subroutine delloc_real1(a)
    implicit none
    real(dp), allocatable, dimension(:), intent(inout) :: a

    if(allocated(a)) deallocate(a)
end subroutine delloc_real1

subroutine delloc_real2(a)
    implicit none
    real(dp), allocatable, dimension(:,:), intent(inout) :: a

    if(allocated(a)) deallocate(a)
end subroutine delloc_real2

subroutine delloc_real3(a)
    implicit none
    real(dp), allocatable, dimension(:,:,:), intent(inout) :: a

    if(allocated(a)) deallocate(a)
end subroutine delloc_real3


subroutine delloc_complex1(a)
    implicit none
    complex(dp), allocatable, dimension(:), intent(inout) :: a

    if(allocated(a)) deallocate(a)
end subroutine delloc_complex1

subroutine delloc_complex2(a)
    implicit none
    complex(dp), allocatable, dimension(:,:), intent(inout) :: a
    if(allocated(a)) deallocate(a)
end subroutine delloc_complex2

subroutine delloc_int1(a)
    implicit none
    integer, allocatable, dimension(:), intent(inout) :: a
    if(allocated(a)) deallocate(a)
end subroutine delloc_int1

subroutine delloc_int2(a)
    implicit none
    integer, allocatable, dimension(:,:), intent(inout) :: a
    if(allocated(a)) deallocate(a)
end subroutine delloc_int2

subroutine alloc_real1(a, n)
    implicit none
    integer,                            intent(in)    :: n
    real(dp), allocatable, dimension(:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A real 1D array is already allocated."
    allocate(a(n)); a=zero;
end subroutine alloc_real1

subroutine alloc_real2(a, n, m)
    implicit none
    integer,                              intent(in)    :: n, m
    real(dp), allocatable, dimension(:,:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A real 2D array is already allocated."
    allocate(a(n, m)); a=zero;
end subroutine alloc_real2

subroutine alloc_real3(a, n, m, k)
    implicit none
    integer,                                intent(in)    :: n, m, k
    real(dp), allocatable, dimension(:,:,:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A real 3D array is already allocated."
    allocate(a(n, m, k)); a=zero;
end subroutine alloc_real3

subroutine alloc_complex1(a, n)
    implicit none
    integer,                               intent(in)    :: n
    complex(dp), allocatable, dimension(:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A complex 1D array is already allocated."
    allocate(a(n)); a=dzero;
end subroutine alloc_complex1

subroutine alloc_complex2(a, n, m)
    implicit none
    integer,                                 intent(in)    :: n, m
    complex(dp), allocatable, dimension(:,:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A complex 2D array is already allocated."
    allocate(a(n, m)); a=dzero;
end subroutine alloc_complex2

subroutine alloc_complex3(a, n, m, k)
    implicit none
    integer,                                    intent(in)    :: n, m, k
    complex(dp), allocatable, dimension(:,:,:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A complex 3D array is already allocated."
    allocate(a(n, m, k)); a=dzero;
end subroutine alloc_complex3

subroutine alloc_int1(a, n)
    implicit none
    integer,                            intent(in)    :: n
    integer, allocatable, dimension(:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A real 1D array is already allocated."
    allocate(a(n)); a=0;
end subroutine alloc_int1

subroutine alloc_int2(a, n, m)
    implicit none
    integer,                              intent(in)    :: n, m
    integer, allocatable, dimension(:,:), intent(inout) :: a
    if(allocated(a).eqv.(.TRUE.)) stop "A real 2D array is already allocated."
    allocate(a(n, m)); a=0;
end subroutine alloc_int2

subroutine reverse_array(arr)
    implicit none
    real(dp), intent(inout) :: arr(:)
    real(dp)                :: rr
    integer                :: ir, ir2, n
    n=size(arr)-1
    ir =1
    ir2=n
    do 
        if(ir>=ir2) exit
        rr = arr(ir)
        arr(ir) = arr(ir2)
        arr(ir2) = rr
        ir = ir + 1
        ir2=ir2-1
    end do
end subroutine reverse_array


subroutine interface2layer(ar1, ar2)
    implicit none
    integer :: ir, n1, n2
    real(dp) :: ar1(:), ar2(:)
    n1 = size(ar1)
    n2 = size(ar2) !! special case where the n2  =  n1 + 1

    do ir = 2, n1
        ar2(ir) = 0.5_dp * (ar1(ir-1) + ar1(ir))
    end do
    ar2(1 ) = 2.0_dp * ar1(1 ) - ar2(2 )
    ar2(n2) = 2.0_dp * ar1(n1) - ar2(n1)
end subroutine interface2layer

end module mod_plib