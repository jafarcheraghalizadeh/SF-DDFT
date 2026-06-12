!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    GLOBAL VARIABLES (SYSTEM DIMENSIONS ETC.)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Global
!
  implicit none
  double precision, parameter :: pi = 3.141592653589793d0
!
  integer, parameter          :: monomer_types=2

  integer,parameter           :: gridx=16, gridy=16, gridz=16
    integer,parameter            :: ngrid = gridx*gridy*gridz
  double precision            :: ds 
  double precision            :: sizex, sizey, sizez, volume, dvol
  integer, parameter  	      :: mixing_dim = 2
  double precision    	      ::KT=0.8d0
  double precision,dimension(2):: norm
    double precision:: lambda =0.5d0
  double precision   :: lambda_max=0.8d0,lambda_min=0.1d0
end module Global
