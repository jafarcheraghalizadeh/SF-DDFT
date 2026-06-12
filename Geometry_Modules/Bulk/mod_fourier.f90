module Fourier_fftw3
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS FOR FOURIER TRANSFORMS - BULK GEOMETRY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use Global
  implicit none

  private

  include "fftw3.f"
 
  integer*8           ::  ftplan,ftbackwardplan
  ! for fortran complex values are supported by fftw, !!!
  ! fftw's row major is warped to column major in Fortran!!!

  double precision, dimension(gridx,gridy,gridz)    :: ftin 
  double complex, dimension(gridx/2+1, gridy,gridz) :: ftout
  double precision, dimension(gridx,gridy,gridz)    :: out
  double precision, dimension(gridx,gridy,gridz)    :: ftbackwardout 
  double complex, dimension(gridx/2+1, gridy,gridz) :: ftbackwardin

  double complex, save, dimension(gridx/2+1,gridy, gridz) :: laplace
  double complex, save, dimension(gridx/2+1,gridy, gridz) :: explaplace
  double complex, save, dimension(gridx/2+1,gridy, gridz,3) :: nabla
!
  public :: Initialize_Fourier
  public :: End_Fourier
  public :: Fourier2Real
  public :: Real2Fourier
  public :: Calculate_Volume_Parameters
  public :: Calculate_Laplace

  public :: nabla
  public :: laplace
  public :: explaplace

contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    INITIALIZATION ROUTINE
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Fourier

    call dfftw_plan_dft_r2c_3d(ftplan,gridx,gridy,gridz,ftin,ftout,FFTW_PATIENT)
    call dfftw_plan_dft_c2r_3d(ftbackwardplan,gridx,gridy,gridz,ftbackwardin,ftbackwardout,FFTW_PATIENT)

  end subroutine Initialize_Fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    CLOSING ROUTINE
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine End_Fourier

    call dfftw_destroy_plan(ftplan)
    call dfftw_destroy_plan(ftbackwardplan)

  end subroutine End_Fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    FOURIER TRANSFORM (BACKWARD)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Fourier2Real(prop,propfourier)

    integer                         :: x,y,z
    double precision, dimension(gridx,gridy,gridz), intent(out)    :: prop
    double complex, dimension((gridx/2+1),gridy,gridz), intent(in) :: propfourier
!
    ftbackwardin = propfourier
    call dfftw_execute(ftbackwardplan)
    prop = ftbackwardout
!
    return
  end subroutine Fourier2Real
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    FOURIER TRANSFORM (FORWARD)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Real2Fourier(prop,propfourier)

    integer                         :: x,y,z
    double precision, dimension(gridx,gridy,gridz), intent(in)      :: prop
    double complex, dimension((gridx/2+1),gridy,gridz), intent(out) :: propfourier
!  
    ftin = prop
    call dfftw_execute(ftplan)
    propfourier = ftout/(gridx*gridy*gridz)
!
    return
  end subroutine Real2Fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    ROUTINE TO TEST FOURIER TRANSFORMS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Test_Fourier
!
    use mtmod

    integer             :: x,y,z
    double precision, dimension(gridx,gridy,gridz)   :: testfield
    double precision, dimension(gridx,gridy,gridz)   :: testfieldnew
    double complex, dimension(gridx/2+1,gridy,gridz) :: testfourier

     do z=1,gridz
      do y=1,gridy
        do x=1,gridx
             testfield(x,y,z)=1.0d0*(1.0d0-2.0d0*grnd())
        end do
      end do
     end do
     do z=1,gridz
      do y=1,gridy
        write(10,'(8g12.4)') (testfield(x,y,z),x=1,8)
      end do
     end do
!
    call real2fourier(testfield,testfourier)
    call fourier2real(testfieldnew,testfourier)
!
     do z=1,gridz
      do y=1,gridy
         write(12,'(8g12.4)') (testfieldnew(x,y,z),x=1,8)
      end do
     end do

  end subroutine Test_Fourier
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE: CALCULATE VOLUME PARAMETERS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Volume_Parameters

    volume = sizex*sizey*sizez
    dvol   = volume/dble(gridx*gridy*gridz) 

  end subroutine Calculate_Volume_Parameters
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE: CALCULATE EXP(LAPLACE*DT) IN FOURIER SPACE
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Laplace

    integer :: x,y,z
    double precision :: qx,qy,qz
    double precision :: dx,dy,dz
!
    dx = sizex/dble(gridx)
    dy = sizey/dble(gridy)
    dz = sizez/dble(gridz)
!
     do z=1,gridz !remember indexshift in fortran
       do y=1,gridy !remember indexshift in fortran
          do x=1,(gridx/2+1) !remember indexshift in fortran
             qx = 2.0d0*pi * dble(x-1) / sizex
             qy = 2.0d0*pi * dble(y-1) / sizey
             qz = 2.0d0*pi * dble(z-1) / sizez
             if (y .gt. (gridy/2+1)) then
               qy = 2.0d0*pi * (dble(y-1)-dble(gridy)) / sizey
             end if
             if (z .gt. (gridz/2+1)) then
               qz = 2.0d0*pi * (dble(z-1)-dble(gridz)) / sizez
             end if
             laplace(x,y,z) = dcmplx(- (qx*qx+qy*qy+qz*qz), 0.0d0)
             explaplace(x,y,z) = dcmplx(exp(-ds*(qx*qx+qy*qy+qz*qz)), 0.0d0)
            ! nabla(x,y,z,1) = dcmplx(0.0d0,(qx))
             !nabla(x,y,z,2) = dcmplx(0.0d0,(qy))
             !nabla(x,y,z,3) = dcmplx(0.0d0,(qz))
             nabla(x,y,z,1) = dcmplx(0.0d0,sin(qx*dx)/dx)
             nabla(x,y,z,2) = dcmplx(0.0d0,sin(qy*dy)/dy)
             nabla(x,y,z,3) = dcmplx(0.0d0,sin(qz*dz)/dz)
            ! sin(qx/dx)/dx
          end do
       end do
     end do

    return

  end subroutine Calculate_Laplace
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Fourier_fftw3
