module functions
contains

! ---------------------------------------------------------------------------------

subroutine make_mesh(x)
use def
double precision, intent(out) :: x(Nx+2*NGC,Nx+2*NGC,2)
integer :: i, j

! Initialize the mesh
do i=1,Nx+2*NGC
   do j=1,Nx+2*NGC
      x(i,j,1)=xmin1+dx*(dble(i-NGC)-half)
      x(i,j,2)=xmin1+dx*(dble(j-NGC)-half)
   enddo
enddo

end subroutine make_mesh

! ---------------------------------------------------------------------------------

subroutine initialization(x,w)
  use def
double precision, intent(in)  :: x(Nx+2*NGC,Nx+2*NGC,2)
double precision, intent(out) :: w(Nx+2*NGC,Nx+2*NGC,nw)

! BEWARE : MUST BE IN PRIMITIVE VARIABLES...
if (iniCond=='Alfven') then
   w(:,:,rho_)=one
   w(:,:,v1_) =zero
   w(:,:,v2_) =0.1d0*dsin(two*dpi*x(:,:,1))
   w(:,:,v3_) =0.1d0*dcos(two*dpi*x(:,:,1))
   w(:,:,p_)  =0.1d0
   w(:,:,rad_)=one
elseif (iniCond=='shock') then
   where (x(:,:,1)<0.d0)
      w(:,:,rho_)=one
      w(:,:,v1_) =zero
      w(:,:,v2_) =zero
      w(:,:,v3_) =zero
      w(:,:,p_)  =one
      w(:,:,rad_)=one
   elsewhere
      w(:,:,rho_)=0.125d0
      w(:,:,v1_) =zero
      w(:,:,v2_) =zero
      w(:,:,v3_) =zero
      w(:,:,p_)  =0.1d0
      w(:,:,rad_)=one
   end where
elseif (iniCond=='2D_Alfven') then
   w(:,:,rho_)= one
   w(:,:,v1_) =-0.1d0*dsin(two*dpi*(x(:,:,1)*dcos(alpha)+x(:,:,2)*dsin(alpha)))*dsin(alpha)
   w(:,:,v2_) = 0.1d0*dsin(two*dpi*(x(:,:,1)*dcos(alpha)+x(:,:,2)*dsin(alpha)))*dcos(alpha)
   w(:,:,v3_) = 0.1d0*dcos(two*dpi*(x(:,:,1)*dcos(alpha)+x(:,:,2)*dsin(alpha)))
   w(:,:,p_)  = 0.1d0
   w(:,:,rad_)= one
else
   stop "Unknown initial conditions"
endif

! ... AND NOW WE CONVERT IN CONSERVATIVE VARIABLES
w(:,:,rho_)=w(:,:,rho_)
w(:,:,m1_) =w(:,:,rho_)*w(:,:,v1_)
w(:,:,m2_) =w(:,:,rho_)*w(:,:,v2_)
w(:,:,m3_) =w(:,:,rho_)*w(:,:,v3_)
w(:,:,e_)  =w(:,:,p_)/(gamma-one)+half*(w(:,:,m1_)**two+w(:,:,m2_)**two+w(:,:,m3_)**two)/w(:,:,rho_)
w(:,:,rad_)=w(:,:,rad_)

end subroutine initialization

! ---------------------------------------------------------------------------------
subroutine computeTimestep(x,w,dt)
  !include 'def.f90'
  use def
double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2)
double precision, intent(in) :: w(Nx+2*NGC,Nx+2*NGC,nw)
double precision, intent(out):: dt
double precision             :: cmax(Nx+2*NGC,Nx+2*NGC)

! BEWARE! THE CFL MUST BE DIVIDED BY THE DIMENSIONALITY!
call getcmax(w,1,cmax)
dt=(dx/maxval(cmax))*cfl/2.d0
call getcmax(w,2,cmax)
dt=min(dt,(dx/maxval(cmax))*cfl/2.d0)

end subroutine computeTimestep

! -----------------------------------------------------------------------------------

subroutine getcmax(w,idir,cmax)
  !include 'def.f90'
  use def
double precision, intent(in) :: w(Nx+2*NGC,Nx+2*NGC,nw)
integer, intent(in) :: idir
double precision, intent(out):: cmax(Nx+2*NGC,Nx+2*NGC)
double precision             ::   cP(Nx+2*NGC,Nx+2*NGC), cM(Nx+2*NGC,Nx+2*NGC)
double precision             :: pthermal(Nx+2*NGC,Nx+2*NGC)
integer :: i, j

pthermal(:,:)=(gamma-one)*(w(:,:,e_)-half*(w(:,:,m1_)**two+w(:,:,m2_)**two+w(:,:,m3_)**two)/w(:,:,rho_))

! we need to compute only the fast magnetoacoustic wave
! since it is the fastest
if (idir==1) then
   cP(:,:)=w(:,:,m1_)/w(:,:,rho_)+&
      dsqrt(half*(        (gamma*pthermal(:,:))/w(:,:,rho_) + &
                   dsqrt(((gamma*pthermal(:,:))/w(:,:,rho_))**two ) ) )
elseif (idir==2) then
   cP(:,:)=w(:,:,m2_)/w(:,:,rho_)+&
      dsqrt(half*(        (gamma*pthermal(:,:))/w(:,:,rho_) + &
                   dsqrt(((gamma*pthermal(:,:))/w(:,:,rho_))**two) ) )
endif

do i=1,Nx+2*NGC
   do j=1,Nx+2*NGC
      cmax(i,j)=abs(cP(i,j))
   enddo
enddo

end subroutine getcmax

! -----------------------------------------------------------------------------------

subroutine integrateInTime(x,w,dt)
use def
use fld
double precision, intent(in)    :: x(Nx+2*NGC,Nx+2*NGC,2), dt
double precision, intent(inout) :: w(Nx+2*NGC,Nx+2*NGC,nw)
! double precision :: f_llf(Nx+2*NGC,Nx+2*NGC,nw) ! fluxes
double precision :: f_llf_X(Nx+2*NGC,Nx+2*NGC,nw) ! fluxes
double precision :: f_llf_Y(Nx+2*NGC,Nx+2*NGC,nw) ! fluxes
integer :: idir

call RiemannFlux(w,1,f_llf_X)
call RiemannFlux(w,2,f_llf_Y)
! w(1+NGC:Nx+NGC,1+NGC:Nx+NGC,:)=w(1+NGC:Nx+NGC,1+NGC:Nx+NGC,:)&
!    -(dt/dx)*(f_llf_X(1+NGC:Nx+NGC,1+NGC:Nx+NGC,:)-f_llf_X(NGC:Nx+NGC-1,1+NGC:Nx+NGC,:))&
!    -(dt/dx)*(f_llf_Y(1+NGC:Nx+NGC,1+NGC:Nx+NGC,:)-f_llf_Y(1+NGC:Nx+NGC,NGC:Nx+NGC-1,:))
w(1+NGC:Nx+NGC,1+NGC:Nx+NGC,1:e_)=w(1+NGC:Nx+NGC,1+NGC:Nx+NGC,1:e_)&
   -(dt/dx)*(f_llf_X(1+NGC:Nx+NGC,1+NGC:Nx+NGC,1:e_)-f_llf_X(NGC:Nx+NGC-1,1+NGC:Nx+NGC,1:e_))&
   -(dt/dx)*(f_llf_Y(1+NGC:Nx+NGC,1+NGC:Nx+NGC,1:e_)-f_llf_Y(1+NGC:Nx+NGC,NGC:Nx+NGC-1,1:e_))

call ImplicitStep(x,w,10)

print*, w(5,5,rad_)

call getBC(x,w)

end subroutine integrateInTime

! -----------------------------------------------------------------------------------

subroutine getBC(x,w)
  !include 'def.f90'
  use def
double precision, intent(in)    :: x(Nx+2*NGC,Nx+2*NGC,2)
double precision, intent(inout) :: w(Nx+2*NGC,Nx+2*NGC,nw)
double precision                :: w0(Nx+2*NGC,Nx+2*NGC,nw)

if (BC=='periodic') then
   ! Periodic BC
   w(1:NGC,:,:)            =w(Nx+1:Nx+NGC,:,:)
   w(Nx+NGC+1:Nx+2*NGC,:,:)=w(NGC+1:NGC+NGC,:,:)
   w(:,1:NGC,:)            =w(:,Nx+1:Nx+NGC,:)
   w(:,Nx+NGC+1:Nx+2*NGC,:)=w(:,NGC+1:NGC+NGC,:)
elseif (BC=='ini') then
   call initialization(x,w0)
   w(1:NGC,:,:)            =w0(1:NGC,:,:)
   w(Nx+NGC+1:Nx+2*NGC,:,:)=w0(Nx+NGC+1:Nx+2*NGC,:,:)
   w(:,1:NGC,:)            =w0(:,1:NGC,:)
   w(:,Nx+NGC+1:Nx+2*NGC,:)=w0(:,Nx+NGC+1:Nx+2*NGC,:)
else
   stop "BC unknown"
endif

end subroutine getBC

! -----------------------------------------------------------------------------------

subroutine RiemannFlux(w,idir,f_llf)
  !include 'def.f90'
  use def
double precision, intent(inout) :: w(Nx+2*NGC,Nx+2*NGC,nw)
integer, intent(in)             :: idir
double precision, intent(out)   :: f_llf(Nx+2*NGC,Nx+2*NGC,nw) ! fluxes @ right edge
double precision                :: f(Nx+2*NGC,Nx+2*NGC,nw)
double precision                :: cmax(Nx+2*NGC,Nx+2*NGC)
integer :: i, j

call getAdvectiveFlux(w,idir,f)

! do idir=1,ndim
! idir=1
call getcmax(w,idir,cmax)
! enddo

if (idir==1) then
   do i=1,Nx+2*NGC-1
      do j=1,Nx+2*NGC-1
         cmax(i,j)=max(cmax(i,j),cmax(i+1,j))
      enddo
   enddo
   do i=1,nw
      f_llf(NGC:Nx+NGC,NGC:Nx+NGC,i)=half*(f(1+NGC:1+Nx+NGC,NGC:Nx+NGC,i)+f(NGC:Nx+NGC,NGC:Nx+NGC,i))&
         -half*cmax(NGC:Nx+NGC,NGC:Nx+NGC)*(w(1+NGC:1+Nx+NGC,NGC:Nx+NGC,i)-w(NGC:Nx+NGC,NGC:Nx+NGC,i))
   enddo
elseif (idir==2) then
   do j=1,Nx+2*NGC-1
      do i=1,Nx+2*NGC-1
         cmax(i,j)=max(cmax(i,j),cmax(i,j+1))
      enddo
   enddo
   do i=1,nw
      f_llf(NGC:Nx+NGC,NGC:Nx+NGC,i)=half*(f(NGC:Nx+NGC,1+NGC:1+Nx+NGC,i)+f(NGC:Nx+NGC,NGC:Nx+NGC,i))&
         -half*cmax(NGC:Nx+NGC,NGC:Nx+NGC)*(w(NGC:Nx+NGC,1+NGC:1+Nx+NGC,i)-w(NGC:Nx+NGC,NGC:Nx+NGC,i))
   enddo
endif

end subroutine RiemannFlux

! -----------------------------------------------------------------------------------

subroutine getAdvectiveFlux(w,idir,f)
  !include 'def.f90'
  use def
double precision, intent(in) :: w(Nx+2*NGC,Nx+2*NGC,nw)
integer, intent(in)          :: idir
double precision, intent(out):: f(Nx+2*NGC,Nx+2*NGC,nw) ! fluxes @ cell centers
double precision             :: wprim(Nx+2*NGC,Nx+2*NGC,nw)

call primitive(w,wprim)

if (idir==1) then
   f(:,:,rho_)=wprim(:,:,rho_)*wprim(:,:,v1_)
   f(:,:,m1_ )=wprim(:,:,rho_)*wprim(:,:,v1_)**two+wprim(:,:,p_)
   f(:,:,m2_ )=wprim(:,:,rho_)*wprim(:,:,v1_)*wprim(:,:,v2_)
   f(:,:,m3_ )=wprim(:,:,rho_)*wprim(:,:,v1_)*wprim(:,:,v3_)
   f(:,:,e_  )=wprim(:,:,v1_)*(w(:,:,e_)+wprim(:,:,p_))
elseif(idir==2) then
   f(:,:,rho_)=wprim(:,:,rho_)*wprim(:,:,v2_)
   f(:,:,m1_ )=wprim(:,:,rho_)*wprim(:,:,v2_)*wprim(:,:,v1_)
   f(:,:,m2_ )=wprim(:,:,rho_)*wprim(:,:,v2_)**two+wprim(:,:,p_)
   f(:,:,m3_ )=wprim(:,:,rho_)*wprim(:,:,v2_)*wprim(:,:,v3_)
   f(:,:,e_  )=wprim(:,:,v2_)*(w(:,:,e_)+wprim(:,:,p_)+half)
endif

end subroutine getAdvectiveFlux

! -----------------------------------------------------------------------------------

subroutine primitive(w,wprim)
  !include 'def.f90'
  use def
double precision, intent(in) :: w(Nx+2*NGC,Nx+2*NGC,nw)
double precision, intent(out):: wprim(Nx+2*NGC,Nx+2*NGC,nw)

wprim(:,:,rho_)=w(:,:,rho_)
wprim(:,:,v1_ )=w(:,:,m1_)/w(:,:,rho_)
wprim(:,:,v2_ )=w(:,:,m2_)/w(:,:,rho_)
wprim(:,:,v3_ )=w(:,:,m3_)/w(:,:,rho_)
wprim(:,:,p_  )=(gamma-one)*(w(:,:,e_)-half*(w(:,:,m1_)**two+w(:,:,m2_)**two+w(:,:,m3_)**two)/w(:,:,rho_))
wprim(:,:,rad_)=w(:,:,rad_)

end subroutine primitive

! ---------------------------------------------------------------------------------

subroutine CV_and_conquer(x,w)
!include 'def.f90'
use def
double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2)
double precision, intent(out) :: w(Nx+2*NGC,Nx+2*NGC,nw)
double precision :: w0(Nx+2*NGC,Nx+2*NGC,nw)
double precision :: Linf, EOC, Linf_old
logical :: file_exists

! Converge (EOC) and accuracy (Linf) measurements
if (.not. abs(tmax-nint(tmax))<smalldble) then
   stop 'Stop! The accuracy analysis requires to know the exact solution, which is the initial state only for integer times'
endif

call initialization(x,w0)

Linf=maxval(abs(w0-w))
inquire(file=EOC_file, exist=file_exists)
if (file_exists) then
   open(1, file=EOC_file, status="old", position="append")
   backspace(1) ! rewind one line
   read(1,*) EOC , EOC, EOC, Linf_old, EOC
   EOC=dlog(Linf/Linf_old)/dlog(dble(Nx/two)/dble(Nx))
else
   open(1, file=EOC_file, status="new", action="write")
   write(1,*) 'AMRlvl N dx Linf EOC'
end if
write(1,*) AMRlvl, Nx, dx, Linf, EOC
close(1)

end subroutine CV_and_conquer

! -----------------------------------------------------------------------------------

! Computes the time ellapsed from start
subroutine chrono(start,start_mess)
  !include 'def.f90'
  use def
double precision, intent(in) :: start
character(len=8), intent(in) :: start_mess ! to have the name of the chrono
double precision :: finish
character(len=200) :: duration

call cpu_time(finish)
write(duration,'(F7.1)') finish-start
duration='It took '//trim(duration)//' sec since '//trim(start_mess)
call followup(duration)

end subroutine chrono

! -----------------------------------------------------------------------------------

subroutine followup(message)
  !include 'def.f90'
  use def
character(LEN=*), intent(in) :: message
! open(4,file=log_file,position="append")
! write(4,*), message
! close(4)
print*, message

end subroutine followup

! -----------------------------------------------------------------------------------

end module functions
