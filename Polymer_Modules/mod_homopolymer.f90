!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR HOMOPOLYMERS 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Molecule_Homopolymer
!
  use Global
  use Propagate
!
  implicit none
!
  private
!
  type homopolymer
    integer          :: label
    integer          :: monomer_type
    double precision :: length
    integer          :: ensemble       ! canonical:1, grand canonical:2
    double precision :: rho            ! total polymer concentration 
    double precision :: mu             ! chemical potential of polymer
    double precision :: energy         ! free energy in SCF field
    integer          :: segments
    double precision :: Qsum
    double precision, dimension(gridx,gridy,gridz)    :: distribution
    double precision, dimension(gridx,gridy,gridz)    :: density
  end type homopolymer
!
  public :: homopolymer
  public :: Initialize_Homopolymer
  public :: Propagate_Homopolymer
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZES HOMOPOLYMER
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Homopolymer(polymer,label,monomer_type,length, &
                                    ensemble,rho,mu)
!
   integer, intent(in)          :: label, ensemble, monomer_type
   double precision, intent(in) :: length, rho, mu
!
   type(homopolymer), intent(out):: polymer
!
   character*15                  :: string_ensemble
!
   if (ensemble.eq.1) then
      string_ensemble = 'canonical'
   else if (ensemble.eq.2) then
      string_ensemble = 'grand canonical'
   else
      print*, 'Invalid Ensemble'
   end if
!
   polymer%label        = label
   polymer%length       = length
   polymer%monomer_type = monomer_type
   polymer%ensemble     = ensemble
   polymer%rho          = rho
   polymer%mu           = mu

   polymer%segments     = nint(polymer%length/ds)
!
   print*, 'Homopolymer ', polymer%label, ' is initialized (', &
        &string_ensemble, ' ensemble).'
!
! ----- Initializations for Testing Purposes ---
!
   polymer%density = polymer%rho * polymer%segments*ds &
     * volume / dvol / dble(gridx*gridy*gridz) ! accounts for boundary effects
!
  end subroutine Initialize_Homopolymer
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE PROPAGATOR FOR A HOMOPOLYMER CHAIN OF LENGTH segments*ds
!		DEDUCE SINGLE CHAIN PARTITION FUNCTION Qsum
!               DEDUCE UNNORMALIZED MONOMER DISTRIBUTION
!               DEDUCE MONOMER DENSITY DUE TO THIS CHAIN
!               DEDUCE CHAIN ENERGY OF THIS CHAIN
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Propagate_Homopolymer(polymer,expfield)
!
    type(homopolymer), intent(inout)  :: polymer
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: expfield

    double precision, dimension(:,:,:,:), allocatable :: qq
    double precision, dimension(gridx,gridy,gridz)    :: prop
!
    integer :: iAllocStatus
    integer :: s,mono,monomer_type
!
    allocate(qq(gridx,gridy,gridz,0:polymer%segments),stat=iAllocStatus)
!
    s = 0
    prop=1.0d0
    qq(:,:,:,s)=prop
!
    do s=1,polymer%segments
       monomer_type = polymer%monomer_type
       call Propagate_Step(prop,expfield(:,:,:,monomer_type))
       qq(:,:,:,s)=prop
    end do
!
    polymer%Qsum = dvol * sum(qq(:,:,:,polymer%segments))
!
    mono = polymer%segments
    polymer%distribution=0.5d0*ds*qq(:,:,:,0)*qq(:,:,:,mono)
    do s=1,mono-1
       polymer%distribution = polymer%distribution &
          + ds*qq(:,:,:,s)*qq(:,:,:,mono-s)
    end do
    polymer%distribution = polymer%distribution &
          + 0.5d0*ds*qq(:,:,:,mono)*qq(:,:,:,0)
!
    if (polymer%ensemble.eq.1) then
      polymer%density = polymer%rho*volume * polymer%distribution/polymer%Qsum
      polymer%mu = log(polymer%rho*volume/polymer%Qsum)
      polymer%energy = polymer%rho*volume * (polymer%mu - 1.0d0)
    else if (polymer%ensemble.eq.2) then
      polymer%density = exp(polymer%mu) * polymer%distribution
      polymer%rho = sum(polymer%density)*dvol / volume / (ds*polymer%segments)
      polymer%energy = - polymer%rho*volume
    else
      print*, 'Invalid ensemble; Cannot calculate Density'
      stop
    end if
!
    deallocate(qq)
!
    return

  end subroutine Propagate_Homopolymer
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Molecule_Homopolymer
