ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c							c
c	Finite Difference Code 			  	c
c	Advection B-GRID MPDATA
c	Shallow Water					c 
c	Time discretization:   Forward
c
c       (7/14/14)
c       beta plane 
c	DIRICHLET BOUNDARIES
c	RECTANGULAR BASIN
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program mpdata

	implicit double precision (a-h,o-z)
	parameter (nsx = 160)
	parameter (nsy = 110)

	dimension vold(nsx,nsy),uold(nsx,nsy)
	dimension unew(nsx,nsy),vnew(nsx,nsy)
	dimension utemp(nsx,nsy),vtemp(nsx,nsy)
	dimension v(nsx,nsy),u(nsx,nsy)
	dimension uxx(nsx,nsy),vyy(nsx,nsy)
	dimension vxx(nsx,nsy),uyy(nsx,nsy)

	dimension f1(nsx,nsy),f2(nsx,nsy)
	dimension hold(0:nsx,0:nsy),htemp(0:nsx,0:nsy)
	dimension hnew(nsx,nsy)

	dimension fau(nsy), fr(nsy)
	real*4 udump(nsx,nsy),vdump(nsx,nsy),hdump(nsx,nsy)

      data ep / .000000001d0 / 

      DONOR(y1, y2, a) = dmax1(0.0d0, a) * y1 + dmin1(0.0d0, a) * y2
      vdyf(x1, x2, a) = (dabs(a) - a*a ) * (dabs(x2) - dabs(x1))
     1      / (dabs(x2) + dabs(x1) + ep)

      vcorr(a, b, y1, y2) = -0.125 * a * b * y1 / (y2)

	open(unit=10,file='mpp_in', status='old')
	read(10,*) nx,ny
	read(10,*) rlx,rly
	read(10,*) 
	read(10,*) amin
	read(10,*) year,day,qtime 
	read(10,*) 
	read(10,*) h0
	read(10,*) 
	read(10,*) feps
	read(10,*) beta
	read(10,*) alpha 
        read(10,*) 
	read(10,*) alambda,agamma
	read(10,*) grav
	read(10,*)

	rewind(10)

c        open (unit=14,file='rs_log')
        open (unit=15,file='rs_h')

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Initialize constants and fields
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	pi  = 2.d0*dacos(-1.0d0)
	dx  = 1.d0*rlx/nx
	dy  = 1.d0*rly/ny
	d1x = 1.d0/dx
	d1y = 1.d0/dy
	dt  = amin*60.d0
	t   = 0.0
	ng  = 1
	nw  = 0
        nt  = 0
	kr  = 0
	a13 = 1.d0/3.d0
	agamma = 10d9*d1x*d1x*d1x*d1x*agamma
	rcu = 0.5d0*dt/dx

        do j = 1,ny+1
          do i = 1,nx+1
            uold(i,j)  = 0.0
            vold(i,j)  = 0.0
            unew(i,j)  = 0.0
            vnew(i,j)  = 0.0
            utemp(i,j) = 0.0
            vtemp(i,j) = 0.0
            uxx(i,j)   = 0.0
            vxx(i,j)   = 0.0
            uyy(i,j)   = 0.0
            vyy(i,j)   = 0.0
            u(i,j)     = 0.0
            v(i,j)     = 0.0
	    hnew(i,j)  = h0
	    f1(i,j)    = 0.0
	    f2(i,j)    = 0.0
         end do
        end do

        do j = 0,ny+1
          do i = 0,nx+1
	    hold(i,j) = h0
	    htemp(i,j)= h0
         end do
        end do

	  do j = 1,ny
                   ay = rly*(1.0d0*(j-1))/(1.d0*ny)
               fau(j) = -alpha*dcos(pi*ay/rly)
               fr(j)  =  feps*(1.0d0+beta*(ay-.5d0*rly)/rly)
	 end do

	write(*,*) nx,ny
	write(*,*)  "Standard Scheme"
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 5	continue

  	do j=1,ny
	  do i=1,nx
	  u(i,j) = rcu*( 1.5*(uold(i,j)+uold(i,j+1) )
     &      -.5*(utemp(i,j)+utemp(i,j+1))   )
	  v(i,j) = rcu*( 1.5*(vold(i,j)+vold(i+1,j) )
     &      -.5*(vtemp(i,j)+vtemp(i+1,j))   )
	  end do
	 end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Up-Wind
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do j=1,ny
         do i=1,nx
               f1(i,j) = donor(hold(i-1,j), hold(i,j), u(i,j))
               f2(i,j) = donor(hold(i,j-1), hold(i,j), v(i,j))
         end do
        end do

         do j=1,ny
            do i=1,nx
               htemp(i,j) = hold(i,j) 
     &         - (f1(i+1,j) - f1(i,j) + f2(i,j+1) - f2(i,j))
         end do
        end do


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Correction Step
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do i=1,nx
	  htemp(i,0)    = htemp(i,1)
	  htemp(i,ny+1) = htemp(i,ny)
	end do

	do j=1,ny
	  htemp(0,j)    = htemp(1,j)
	  htemp(nx+1,j) = htemp(nx,j)
	end do

         do j=1,ny
            do  i=1,nx
              u(i, j) = vdyf( htemp(i-1,j),htemp(i,j),u(i, j) )
     &        + vcorr( u(i,j),4.0*rcu*( 1.5*(vold(i,j)+vold(i,j+1) )
     &                                -.5*(vtemp(i,j)+vtemp(i,j+1))),
     &        dabs(htemp(i-1, j+1)) + dabs(htemp(i, j+1)) -
     &        dabs(htemp(i-1, j-1)) - dabs(htemp(i, j-1)),
     &        dabs(htemp(i-1, j+1)) + dabs(htemp(i, j+1)) +
     &        dabs(htemp(i-1, j-1)) + dabs(htemp(i, j-1)) + ep)

              v(i, j) = vdyf( htemp(i,j-1),htemp(i,j),v(i,j) )
     &          + vcorr( v(i,j),4.0*rcu*( 1.5*(uold(i,j)+uold(i+1,j) )
     &                                  -.5*(utemp(i,j)+utemp(i+1,j))),
     &          dabs(htemp(i+1, j-1)) + dabs(htemp(i+1, j)) -
     &          dabs(htemp(i-1, j-1)) - dabs(htemp(i-1, j)),
     &          dabs(htemp(i+1, j-1)) + dabs(htemp(i+1, j)) +
     &          dabs(htemp(i-1, j-1)) + dabs(htemp(i-1, j)) + ep )
         end do
        end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Up-Wind AGAIN
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do j=1,ny
         do i=1,nx
               f1(i, j) = donor(htemp(i-1, j), htemp(i, j), u(i, j))
               f2(i, j) = donor(htemp(i, j-1), htemp(i, j), v(i, j))
         end do
        end do

         do j=1,ny
            do i =1,nx
               hnew(i,j) = htemp(i,j) 
     &         - (f1(i+1,j) - f1(i,j) + f2(i,j+1) - f2(i,j))
         end do
        end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Up-date Velocities
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do j=1,ny+1
	   do i=1,nx+1
	    u(i,j) = 1.5*uold(i,j)-.5*utemp(i,j)
	    v(i,j) = 1.5*vold(i,j)-.5*vtemp(i,j)
	   end do
	end do

        do j=2,ny
         do i=2,nx

	    uux = u(i,j)*(u(i+1,j)-u(i-1,j))
	    vuy = v(i,j)*(u(i,j+1)-u(i,j-1))
	    uvx = u(i,j)*(v(i+1,j)-v(i-1,j))
	    vvy = v(i,j)*(v(i,j+1)-v(i,j-1))

	    uxx(i,j) = u(i+1,j)-2.0*u(i,j)+u(i-1,j)
	    vxx(i,j) = v(i+1,j)-2.0*v(i,j)+v(i-1,j)
	    uyy(i,j) = u(i,j+1)-2.0*u(i,j)+u(i,j-1)
	    vyy(i,j) = v(i,j+1)-2.0*v(i,j)+v(i,j-1)

	  rpx = .5*(hnew(i,j-1)+hnew(i,j)-hnew(i-1,j-1)-hnew(i-1,j))
	  rpy = .5*(hnew(i-1,j)+hnew(i,j)-hnew(i-1,j-1)-hnew(i,j-1))

	   su = 1.5*uold(i,j)-.5*utemp(i,j)
	   sv = 1.5*vold(i,j)-.5*vtemp(i,j)

	  res = -0.5*d1x*uux-0.5*d1y*vuy +fr(j)*sv
     &	      + alambda*(d1x*d1x*uxx(i,j)+d1y*d1y*uyy(i,j)) + fau(j)
     &	      - agamma*(uxxxx+uyyyy)

	 resu = uold(i,j)+dt*(res-d1x*grav*rpx)

	  res = -0.5*d1x*uvx-0.5*d1y*vvy - fr(j)*su
     &	      + alambda*(d1x*d1x*vxx(i,j)+d1y*d1y*vyy(i,j)) 
     &	      - agamma*(vxxxx+vyyyy)

         resv = vold(i,j)+dt*(res-d1y*grav*rpy)

	unew(i,j) = resu
	vnew(i,j) = resv

	 end do
	end do

	t = t + dt

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Writting Statments
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	nw = nw+1
        if(nw*amin.ge.year*525600.0) then
            ng = ng+1
	  if(ng*amin.ge.day*24.0*60.0) then
          do i = 1,nx
           do j = 1,ny
            write(15,*)i, j, hold(i,j)-h0
           end do
          end do
          ng = 0
          nt = nt +1
          write(*,*) nt
	end if
  

	end if




c 666	format(f12.1,3f12.5,i5)
c 667	format(f7.0,2f8.2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Stop or return to loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if (t.ge.qtime*31536000.0) then
	  go to 10
	end if

        do j=1,ny+1
         do i=1,nx+1
	    utemp(i,j) = uold(i,j)
	    vtemp(i,j) = vold(i,j)
	    uold(i,j)  = unew(i,j)
	    vold(i,j)  = vnew(i,j)
	    udump(i,j) = unew(i,j)
	    vdump(i,j) = vnew(i,j)
         end do
        end do

        do j=1,ny
         do i=1,nx
	  hold(i,j)  = hnew(i,j)
	  hdump(i,j) = hnew(i,j)
         end do
        end do

	go to 5
 10	continue

	stop
	end
