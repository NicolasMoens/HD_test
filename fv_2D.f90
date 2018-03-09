program fv_2D
use functions
use def

double precision :: x(Nx+2*NGC,Nx+2*NGC,2)
double precision :: w(Nx+2*NGC,Nx+2*NGC,nw), w0(Nx+2*NGC,Nx+2*NGC,nw)
double precision :: dt, t=0.d0
double precision :: chrono_0
character(len=400) :: out_file
integer :: i, j, iw, itsave=0

! Compilation
! gfortran functions.f90 fv_1D.f90 -o fv_1D

call cpu_time(chrono_0)

call make_mesh(x)

print*, 'Initialization - - - - - - - - - -'

call initialization(x,w)

print*, 'Core mod starts - - - - - - - - - -'

write (out_file, "(A5,I0.4)") "wdat_", itsave
open(1,file=out_file)
   call primitive(w,w0)
   write(1,'(1pe20.12)') t
   do j=1,Nx+2*NGC
      do i=1,Nx+2*NGC
         write(1,'(200(1pe20.12))') x(i,j,1), x(i,j,2), (w0(i,j,iw),iw=1,nw)
      enddo
   enddo
close(1)
itsave=itsave+1

do while (t<tmax)

   call computeTimestep(x,w,dt)
   ! print*, dt, t, tmax

   ! make sure we do not overshoot tmax
   dt=min(dt,tmax-t)

   call integrateInTime(x,w,dt)

   t=t+dt

   ! Save condition
   if (floor(t/dtsave)/=floor((t-dt)/dtsave)) then
      write (out_file, "(A5,I0.4)") "wdat_", itsave
      open(1,file=out_file)
         call primitive(w,w0)
         write(1,'(1pe20.12)') t
         do j=1,Nx+2*NGC
            do i=1,Nx+2*NGC
               write(1,'(200(1pe20.12))') x(i,j,1), x(i,j,2), (w0(i,j,iw),iw=1,nw)
            enddo
         enddo
      close(1)
      itsave=itsave+1
   endif

   print*, 100.d0*(t/tmax),'%', dt

enddo

write (out_file, "(A5,I0.4)") "wdat_", itsave
open(1,file=out_file)
   call primitive(w,w0)
   write(1,'(1pe20.12)') t
   do j=1,Nx+2*NGC
      do i=1,Nx+2*NGC
         write(1,'(200(1pe20.12))') x(i,j,1), x(i,j,2), (w0(i,j,iw),iw=1,nw)
      enddo
   enddo
close(1)

if (CV_analysis) then
   call CV_and_conquer(x,w)
endif

print*, 'Terminated - - - - - - - - - - - -'

call chrono(chrono_0,chrono_0_mess)

end program fv_2D
