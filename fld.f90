module fld
contains

! ---------------------------------------------------------------------------------

subroutine radiationflux(x,w,rf,rp)
use def
double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2), w(Nx+2*NGC,Nx+2*NGC,nw)
double precision, intent(out) :: rf(Nx+2*NGC,Nx+2*NGC,2),rp(Nx+2*NGC,Nx+2*NGC)
double precision :: lambda(Nx+2*NGC,Nx+2*NGC)
double precision :: fld_r(Nx+2*NGC,Nx+2*NGC)
double precision :: gradE(Nx+2*NGC,Nx+2*NGC,2)

!calculate gradient of E
call grad(x(:,:,:), w(:,:,rad_),gradE(:,:,:))

!calculating R = norm(gradE)/(chi*E)
fld_r(:,:) = dsqrt(gradE(:,:,1)**two + gradE(:,:,2)**two)/(kappa*w(:,:,rho_)*w(:,:,rad_))

!calculating Lambda = (2+R)/(6+3R+R**2)
lambda(:,:) = (2 + fld_r(:,:))/(6 + 3*fld_r(:,:) + fld_r(:,:)**two)

!radiationflux given by F_i = c*Lambda/kappa*rho * gradE
rf(:,:,1) = (clight*lambda(:,:))/(kappa*w(:,:,rho_)) * gradE(:,:,1)
rf(:,:,2) = (clight*lambda(:,:))/(kappa*w(:,:,rho_)) * gradE(:,:,2)

! !radiation pressure by P = (l + l**2 R**2) E
! rp(:,:) = (lambda(:,:) + lambda(:,:)**two * fld_r(:,:)**two ) * w(:,:,rad_)

end subroutine radiationflux

! -----------------------------------------------------------------------------------

subroutine grad(x,q,gradq)
use def
double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2), q(Nx+2*NGC,Nx+2*NGC)
double precision, intent(out) :: gradq(Nx+2*NGC,Nx+2*NGC,2)
integer :: i,j

do i=NGC,Nx+NGC
  do j=NGC,Nx+NGC
     gradq(i,j,1) = (q(i+1,j)-q(i-1,j))/(x(i+1,j,1)-x(i-1,j,1))
     gradq(i,j,2) = (q(i,j+1)-q(i,j-1))/(x(i,j+1,2)-x(i,j-1,2))
  enddo
enddo

end subroutine grad

! -----------------------------------------------------------------------------------

subroutine make_matrix(x,w,dw,sweepdir,diag,sub,sup)
use def
integer, intent(in) :: sweepdir
double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2), w(Nx+2*NGC,Nx+2*NGC,nw), dw
double precision, intent(out) :: diag(Nx+2*NGC,Nx+2*NGC),sub(Nx+2*NGC,Nx+2*NGC),sup(Nx+2*NGC,Nx+2*NGC)
double precision :: h
double precision :: beta(Nx+2*NGC)
double precision :: lambda(Nx+2*NGC,Nx+2*NGC)
double precision :: fld_r(Nx+2*NGC,Nx+2*NGC)
double precision :: gradE(Nx+2*NGC,Nx+2*NGC,2)
double precision :: D_center(Nx+2*NGC,Nx+2*NGC)
double precision :: D(Nx+2*NGC,Nx+2*NGC,2)

integer :: i,j

!calculate gradient of E
call grad(x(:,:,:), w(:,:,rad_),gradE(:,:,:))

!calculating R = norm(gradE)/(chi*E)
fld_r(:,:) = dsqrt(gradE(:,:,1)**two + gradE(:,:,2)**two)/(kappa*w(:,:,rho_)*w(:,:,rad_))

!calculating Lambda = (2+R)/(6+3R+R**2)
lambda(:,:) = (2 + fld_r(:,:))/(6 + 3*fld_r(:,:) + fld_r(:,:)**two)

!calculate diffusion coefficient
D_center(:,:) = clight*lambda(:,:)/(kappa*w(:,:,rho_))

!> Go from cell center to cell face
do i = 2, NX+2*NGC
do j = 2, NX+2*NGC
  D(i,j,1) = (D_center(i,j) + D_center(i-1,j))/2.d0
  D(i,j,2) = (D_center(i,j) + D_center(i,j-1))/2.d0
enddo
enddo

!calculate h
h = dw/(two*dx**two)

if (sweepdir == 1) then
  !calculate matrix for sweeping in 1-direction
  do j = NGC,Nx+NGC
   !calculate beta
   do i = NGC,Nx+NGC
     beta(i) = one + dw/(two*dt) + h*(D(i+1,j,1))
   enddo

   do i = NGC,Nx+NGC
     diag(i,j) = beta(i)
     sub(i,j) = -h*D(i,j,1)
     sub(i,j) = -h*D(i+1,j,1)
   enddo
   !> Fix boundary conditions on matrix nd stuf
  enddo
elseif ( sweepdir == 2 ) then
  !calculate matrix for sweeping in 2-direction
  do j = NGC,Nx+NGC
   !calculate beta
   do i = NGC,Nx+NGC
     beta(i) = one + dw/(two*dt) + h*(D(j,i+1,2))
   enddo

   do i = NGC,Nx+NGC
     diag(i,j) = beta(i)
     sub(i,j) = -h*D(i,j,2)
     sub(i,j) = -h*D(j,i+1,2)
   enddo
   !> Fix boundary conditions on matrix nd stuf
  enddo
else
  stop "sweepdir unknown"
endif

end subroutine make_matrix

! -----------------------------------------------------------------------------------

subroutine make_bvec(x,w,dw,E_n,E_m,sweepdir,bvec)
use def
 integer, intent(in) :: sweepdir
 double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2), w(Nx+2*NGC,Nx+2*NGC,nw), dw
  double precision, intent(in) :: E_n(Nx+2*NGC,Nx+2*NGC), E_m(Nx+2*NGC,Nx+2*NGC)
 double precision, intent(out) :: bvec(Nx+2*NGC,Nx+2*NGC)
 double precision :: lambda(Nx+2*NGC,Nx+2*NGC)
 double precision :: fld_r(Nx+2*NGC,Nx+2*NGC)
 double precision :: gradE(Nx+2*NGC,Nx+2*NGC,2)
 double precision :: D_center(Nx+2*NGC,Nx+2*NGC)
 double precision :: D(Nx+2*NGC,Nx+2*NGC,2)
 double precision :: h
 integer :: i,j


 !calculate gradient of E
 call grad(x(:,:,:), w(:,:,rad_),gradE(:,:,:))

 !calculating R = norm(gradE)/(chi*E)
 fld_r(:,:) = dsqrt(gradE(:,:,1)**two + gradE(:,:,2)**two)/(kappa*w(:,:,rho_)*w(:,:,rad_))

 !calculating Lambda = (2+R)/(6+3R+R**2)
 lambda(:,:) = (2 + fld_r(:,:))/(6 + 3*fld_r(:,:) + fld_r(:,:)**two)

 !calculate diffusion coefficient
 D_center(:,:) = clight*lambda(:,:)/(kappa*w(:,:,rho_))

 !> Go from cell center to cell face
 do i = 2, NX+2*NGC
 do j = 2, NX+2*NGC
   D(i,j,1) = (D_center(i,j) + D_center(i-1,j))/2.d0
   D(i,j,2) = (D_center(i,j) + D_center(i,j-1))/2.d0
 enddo
 enddo
 !calculate h
 h = dw/(two*dx**two)


if (sweepdir == 1) then
 do j = NGC,Nx+NGC
   do i = NGC,Nx+NGC
     bvec(i,j) = (1 + h*(D(i,j+1,2)+D(i,j,2)))*E_m(i,j) &
     + h*D(i,j+1,2)*E_m(i,j+1) + h*D(i,j,2)*E_m(i,j-1) + dw/(2*dt)*E_n(i,j)
   enddo
 enddo
elseif ( sweepdir == 2 ) then
  do j = NGC,Nx+NGC
    do i = NGC,Nx+NGC
      bvec(i,j) = (1 + h*(D(j+1,i,1)+D(i,j,1)))*E_m(i,j) &
      + h*D(j+1,1,1)*E_m(j+1,1) + h*D(i,j,1)*E_m(j-1,i) + dw/(2*dt)*E_n(i,j)
    enddo
  enddo
else
  stop "sweepdir unknown"
endif
end subroutine make_bvec

! -----------------------------------------------------------------------------------

subroutine solve_tridiag(diag,sub,sup,bvec,E_m)
  use def
  implicit none
  !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
  !	 b - the main diagonal
  !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
  !	 d - right part
  !	 x - the answer
  !	 n - number of equations

  double precision, intent(in) :: diag(Nx+2*NGC), bvec(Nx+2*NGC)
  double precision, intent(in) :: sub(Nx+2*NGC), sup(Nx+2*NGC)
  double precision, intent(out) :: E_m(Nx+2*NGC)
  double precision :: cp(Nx+2*NGC), dp(Nx+2*NGC)
  double precision :: m
  integer :: i
  ! initialize c-prime and d-prime
  cp(1) = sup(1)/diag(1)
  dp(1) = bvec(1)/diag(1)
  ! solve for vectors c-prime and d-prime
  do i = 2,Nx+2*NGC
    m = diag(i)-cp(i-1)*sub(i)
    cp(i) = sup(i)/m
    dp(i) = (bvec(i)-dp(i-1)*sub(i))/m
  enddo
  ! initialize x
  E_m(Nx+2*NGC) = dp(Nx+2*NGC)
  ! solve for x from the vectors c-prime and d-prime
  do i = Nx+2*NGC-1, 1, -1
    E_m(i) = dp(i)-cp(i)*E_m(i+1)
  end do

end subroutine solve_tridiag

! -----------------------------------------------------------------------------------

subroutine ImplicitStep(x,w,w_max)
use def
integer, intent(in) :: w_max
double precision, intent(in) :: x(Nx+2*NGC,Nx+2*NGC,2)
double precision, intent(inout) :: w(Nx+2*NGC,Nx+2*NGC,nw)
double precision :: E_m(Nx+2*NGC,Nx+2*NGC), E_n(Nx+2*NGC,Nx+2*NGC)
double precision :: diag(Nx+2*NGC,Nx+2*NGC), bvec(Nx+2*NGC,Nx+2*NGC)
double precision :: sub(Nx+2*NGC,Nx+2*NGC), sup(Nx+2*NGC,Nx+2*NGC)
double precision :: dw
integer m

E_n(:,:) = w(:,:,rad_)
E_m(:,:) = w(:,:,rad_)

do m = 1,w_max
  !> Set pseudotimestep
  dw = dx/4.d0*((xmax1-xmin1)/dx)**((m-1)/(w_max-1))

  !> Setup matrix and vector for sweeping in direction 1
  call make_matrix(x,w,dw,1,diag,sub,sup)
  call make_bvec(x,w,dw,E_n,E_m,1,bvec)
  do j = NGC,Nx+NGC
    call solve_tridiag(diag(:,j),sub(:,j),sup(:,j),bvec(:,j),E_m(:,j))
  enddo
  !> Setup matrix and vector for sweeping in direction 2
  call make_matrix(x,w,dw,2,diag,sub,sup)
  call make_bvec(x,w,dw,E_n,E_m,2,bvec)
  do j = NGC,Nx+NGC
    call solve_tridiag(diag(:,j),sub(:,j),sup(:,j),bvec(:,j),E_m(:,j))
  enddo
enddo

w(:,:,rad_) = E_m(:,:)

end subroutine ImplicitStep

! -----------------------------------------------------------------------------------

end module fld
