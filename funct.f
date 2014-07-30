c-------------------------------------------------------------------------
c integrate y' = -y+f,  y(0)=y0 to time T
c 6/23/2014
c--------------------------------------------------------------------------

             implicit double precision (a-h,o-z)
             parameter (nt = 250)

             dimension y(nt),f(nt)

             open(unit=1,file='plot.dat')

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Initialize constants and fields
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

             y0 = 2.d0
             dt = .01d0
             T  = 1.d0
             nq = T/dt

        if(nq.ge.nt) then
           write(*,*) 'nt must be larger than'
           write(*,*) nq
           stop
        end if

           do i=1,nq
             f(i) = .1d0
             y(i) = 0.d0
           end do
             y(1) = y0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Forward Euler
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=2,nq
           y(i) = y(i-1)+dt*(f(i-1)-y(i-1))
         end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write Results
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1,nq
           write(1,10) dble((i-1))*dt,y(i)
         end do

 10      format(f4.3,f7.3)

         stop
         end
