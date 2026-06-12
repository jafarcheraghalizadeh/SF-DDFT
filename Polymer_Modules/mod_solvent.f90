!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR SOLVENT
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Molecule_Solvent
!
  use Global
  use Propagate
!
  implicit none
!
  private
!
  type monomeric_solvent
    integer          :: label
    integer          :: monomer_type
    integer          :: ensemble       ! canonical:1, grand canonical:2
    double precision :: rho            ! total solvent concentration 
    double precision :: mu             ! chemical potential of polymer
    double precision :: energy         ! free energy in SCF field
    double precision :: Qsum
    double precision, dimension(gridx,gridy,gridz)    :: distribution
    double precision, dimension(gridx,gridy,gridz)    :: density
  end type monomeric_solvent
!
  public :: monomeric_solvent
  public :: Initialize_Solvent
  public :: Update_Solvent
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZES SOLVENT
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Solvent(solvent,label,monomer_type,ensemble,rho,mu)
!
   integer, intent(in)          :: label, ensemble, monomer_type
   double precision, intent(in) :: rho,mu
!
   type(monomeric_solvent), intent(out):: solvent
!
   character*15                  :: string_ensemble
!
   integer :: iAllocStatus
!
   if (ensemble.eq.1) then
      string_ensemble = 'canonical'
   else if (ensemble.eq.2) then
      string_ensemble = 'grand canonical'
   else
      print*, 'Invalid Ensemble'
   end if
!
   solvent%label        = label
   solvent%monomer_type = monomer_type
   solvent%ensemble     = ensemble
   solvent%rho          = rho
   solvent%mu           = mu
!
   print*, 'Solvent ', solvent%label, ' is initialized (', &
        &string_ensemble, ' ensemble).'
!
! ----- Initializations for Testing Purposes ---
!
   solvent%density = solvent%rho &
      * volume / dvol / dble(gridx*gridy*gridz) ! accounts for boundary effects
!
  end subroutine Initialize_Solvent
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          CALCULATE UNNORMALIZED SOLVENT DISTRIBUTION
!	        DEDUCE SINGLE PARTICLE PARTITION FUNCTION Qsum
!               DEDUCE SOLVENT DENSITY 
!               DEDUCE ENERGY OF SOLVENT
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Update_Solvent(solvent,field)
!
    type(monomeric_solvent), intent(inout)  :: solvent
    double precision, dimension(gridx,gridy,gridz), intent(in)  :: field
!
    solvent%distribution = exp(-field)
    solvent%Qsum = dvol * sum(solvent%distribution)
!
    if (solvent%ensemble.eq.1) then
      solvent%density = solvent%rho*volume * solvent%distribution/solvent%Qsum
      solvent%mu = log(solvent%rho*volume/solvent%Qsum)
      solvent%energy = solvent%rho*volume * (solvent%mu - 1.0d0)
    else if (solvent%ensemble.eq.2) then
      solvent%density = exp(solvent%mu) * solvent%distribution
      solvent%rho = sum(solvent%density)*dvol / volume 
      solvent%energy = - solvent%rho*volume
    else
      print*, 'Invalid ensemble; Cannot calculate Density'
      stop
    end if
!
    return

  end subroutine Update_Solvent
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Molecule_Solvent
