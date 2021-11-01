program test

        implicit none
        integer     ::   a(10000,10000),b(10000,10000),c(10000,10000)
        real        :: start, endtime
        integer     :: i,j
        a=1
        b=1
        do i=1,10000
           do j=1,10000
           a(i,j)=i+j
           b(i,j)=2*i+j
           enddo
        enddo
        call cpu_time(start)
        do i=1,100
           c=a*b
        enddo
        call cpu_time(endtime)
        print *,'c=a*b used',endtime-start,' s'

        call cpu_time(start)
        do i=1,100
           do j=1,10000
           c(:,j)=b(:,j)*a(:,j)
           enddo
        enddo
        call cpu_time(endtime)
        print *,'c(:,i)=a(:,i)*b(:,i)',endtime-start,' s'

        call cpu_time(start)
        do i=1,100
           do j=1,10000
           c(j,:)=b(j,:)*a(j,:)
           enddo
        enddo
        call cpu_time(endtime)
        print *,'c(i,:)=a(i,:)*b(i,:)',endtime-start,' s'

end program test 

