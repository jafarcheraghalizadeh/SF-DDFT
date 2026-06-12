!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           GENERAL PROPAGATOR SUBROUTINES
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Propagate
!
  use Global
!
  implicit none
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE STEP: CALCULATE EXP(-field*ds*0.5) 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Expfield(field,expfield)
!
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)    :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out)   :: expfield

!
      expfield = exp(-0.5d0*ds*field)
 !print*,maxval(expfield)
    return

  end subroutine Calculate_Expfield
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           PROPAGATE prop BY ONE STEP ds
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Propagate_Step(prop,expfield)
!
    use Fourier_fftw3            
     use, intrinsic :: ieee_arithmetic
    implicit none
!
    double precision, dimension(gridx,gridy,gridz), intent(in)    :: expfield
    double precision, dimension(gridx,gridy,gridz), intent(inout) :: prop
    double complex, dimension(gridx/2+1,gridy,gridz)   :: propfourier
    integer                         :: x,y,z
    double precision::hu
!
!   --- expfield = exp(-field*ds/2.) 
!           is evaluated in Calculate_Expfield (this module)
!   --- explaplace  = exp(-q^2*ds) 
!           is evaluated in Calculate_Laplace (module Fourier_fftw3)
    prop = expfield * prop
if (maxval(prop) > huge(hu)) then
print*, 'prop is infinity1'
print*,minval(prop),maxval(prop)
end if
    call Real2Fourier(prop,propfourier)

    propfourier = explaplace * propfourier

    call Fourier2Real(prop,propfourier)

    prop = expfield * prop 
   return

  end subroutine Propagate_Step

end module Propagate
