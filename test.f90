program test
    implicit none
    double precision    ::  a(5,5)
    integer             ::  i,j

    function test(a)
        implicit none
        double precision    ::  a,test
        integer             ::  i,j
        do i=1,size(a,1)
            do j=1,size(a,2)
            a(i,j)=i+j
            write(*,*) a(i,j)
            enddo
        enddo
        test=a
    end function test
    
end program test