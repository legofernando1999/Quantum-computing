module multiplication
    interface operator (.ext.)
        module procedure ExtMul
    end interface
    
    contains
        function ExtMul(u,v) result(w)
            implicit none 
            integer :: i,j,N
            real*8,dimension(:),intent(in) :: u,v
            real*8,dimension(size(u),size(v)) :: w
            N=size(u)
            w=0.0
            do i=1,N
                do j=1,N
                    w(i,j)=u(i)*v(j)
                end do
            end do
        end function ExtMul

end module multiplication

program Grover
    use multiplication
    implicit none
    integer :: i,j,k,l
    integer,parameter :: n=4
    real*8,parameter :: c=1/sqrt(2.0)
    real*8,dimension(2) :: e0,plus
    real*8,dimension(2,2) :: H
    real*8,dimension(2**n) :: S,W,R
    real*8,dimension(2**n,2**n) :: Uw,Us,Id

    e0=0.0;e0(1)=1.0
    W=0.0;W(7)=1.0
    
    Id=0.0
    do i=1,2**n
        Id(i,i)=1.0
    end do  
    
    H=c;H(2,2)=-c
    plus=matmul(H,e0)
    
    S=0.0
    do i=1,2
        S(i)=plus(i)
    end do

    do l=1,n-1
        k=1
        R=S
        do i=1,2
            do j=1,2**l
                S(k)=plus(i)*R(j)
                k=k+1
            end do
        end do
    end do
    
    Uw=Id-2*(W.ext.W)
    Us=2*(S.ext.S)-Id
    
    open(1,file='Histogram.dat')
    R=S
    do j=1,2**n
        write(1,*) j,R(j)
    end do
    write(1,*)
    write(1,*)
    
    do i=1,nint(sqrt(2.0**n))
        R=matmul(Uw,R)
        do j=1,2**n
            write(1,*) j,R(j)
        end do
        write(1,*)
        write(1,*)
        R=matmul(Us,R)
        do j=1,2**n
            write(1,*) j,abs(R(j))
        end do
        write(1,*)
        write(1,*)
    end do
    
    open(2,file='Histogram.gpl')
    write(2,*) 'set xlabel "Estados"'
    write(2,*) 'set ylabel "Amplitud"'
    write(2,*) 'set yrange [-1.2:1.2]'
    write(2,*) 'set style data histogram'
    write(2,*) 'set style histogram cluster gap 1'
    write(2,*) 'set style fill solid'
    write(2,*) 'set boxwidth 0.9'
    write(2,*) 'set terminal gif animate opt delay 100'
    write(2,*) 'set font "Helvetica,30"'
    write(2,*) 'set tics font",',30-n*4,'"'
    write(2,*) 'set output "Histogram.gif"'
    l=nint(sqrt(2.0**(n-1)))*2
    write(2,'(" do for [i=0:",i0,"] {")') l
    write(2,*) '    plot "Histogram.dat" index i using 2:xtic(1) linecolor "blue" notitle'
    write(2,*) '}'

    call execute_command_line('gnuplot Histogram.gpl')
    call execute_command_line('gwenview "Histogram.gif"')
    
    i=1
    do while (i.le.2**n)
        if (nint(R(i)).eq.1) then
            print*, 'El objeto buscado esta en la posicion:',i
            exit
        end if
        i=i+1
    end do
    
end program Grover
