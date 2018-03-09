module def

implicit none
public

double precision, parameter :: smalldble=1.d-10, bigdble=1.d99
double precision, parameter :: one=1.d0, dpi=dacos(-1.d0), zero=0.d0, half=0.5d0, two=2.d0

! MKSA International units system
double precision, parameter :: kboltz = 1.3806488D-23
double precision, parameter :: stefan = 5.670373d-8
double precision, parameter :: mnucl = 1.6733D-27
double precision, parameter :: msun = 1.989D+30
double precision, parameter :: rsolar = 7.0D+8
double precision, parameter :: secinday= 8.64d4
double precision, parameter :: secinyear= 3.1557600d7
double precision, parameter :: ggrav  = 6.67384D-11
double precision, parameter :: clight = 2.99792458d8
double precision, parameter :: Lsun = 3.846d26
double precision, parameter :: AU = 1.49597871d11
double precision, parameter :: Jsun = 1.1d42 ! ang. mom. for 25d
double precision, parameter :: sigmaE = 6.652458734d-29 ! = sigma_thomson_scattering_on_free_electrons
double precision, parameter :: kpE = 3.98d-2!3.1d-2! = sigma_thomson / mass_proton.

double precision, parameter :: kappa = 0.34d0


! For kpE, we used the formula (65) of KUDRITZKI+89 w/ I_HE=2 & N_He/N_H = 0.2

integer, parameter :: AMRlvl=4
integer, parameter :: Nx=10*2**AMRlvl, ndim=1, nw=6, NGC=1 ! NGC is # of ghost cells on each side
double precision, parameter :: xmin1=0.d0, xmax1=dsqrt(two), dx=(xmax1-xmin1)/dble(Nx), cfl=0.45d0, tmax=1.d0, dtsave=tmax/5.d0
double precision, parameter :: gamma=5.d0/3.d0, alpha=dpi/4.d0
logical, parameter :: CV_analysis=.true.
character(len=400), parameter :: BC='periodic', iniCond='2D_Alfven'

integer, parameter :: rho_=1, m1_=rho_+1, m2_=m1_+1, m3_=m2_+1, e_=m3_+1, rad_ =e_+1
integer, parameter :: m0_=rho_, v1_=m1_, v2_=m2_, v3_=m3_, p_=e_

character(len=8), parameter :: chrono_0_mess='chrono_0'

! character(len=400), parameter :: out_file='w.dat'
character(len=400), parameter :: ini_file='w0.dat'
character(len=400), parameter :: EOC_file='EOC.dat'

end module def
