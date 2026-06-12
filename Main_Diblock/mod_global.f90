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
  double precision, parameter  :: pi = 3.141592653589793d0
!
  integer, parameter           :: monomer_types=2
!
  integer,parameter            :: gridx=18, gridy=18, gridz=18
  integer,parameter            :: ngrid = gridx*gridy*gridz
  double precision             :: ds 
  double precision             :: sizex, sizey, sizez, volume, dvol
  double precision,dimension(2):: norm
  double precision             :: lambda=0.40d0
  double precision,parameter   :: lambda_max=1.0,lambda_min=0.1
  integer, parameter           :: mixing_dim=2 ! should at least be 2
  double precision,parameter   :: KT = 0.20d0!  The Noise standard deviation

end module Global
