!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!        ROUTINES FOR FLORY HUGGINS TYPE INTERACTIONS IN SOLUTION
!             IN A MULTICOMPONENT SYSTEM
!             WITH INCOMPRESSIBLE IMPLICIT SOLVENT 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Interactions_FloryHuggins_Solution
!
  use Global
!
  implicit none
!
  private
!
  integer, parameter     :: components = monomer_types     ! solvent is component 0
!
  double precision, save :: kappa               ! Helfand compressibility parameter
!      ---------- ALLOW OVERCOMPRESSING PURE POLYMER MELT FOR STABILITY REASONS
  double precision, save :: monomer_volume(0:components)   ! monomer volumes
  double precision, save :: chi(0:components,0:components) ! interaction parameters
  double precision, save :: chi_eff(components,components) ! interaction parameters
  double precision, save :: zeta                ! rescaled inverse solvent volume
!
  double precision, save :: small_number = 1.0d-3
!
  public :: Initialize_Interactions
  public :: Calculate_Mean_Fields
  public :: Interaction_Energy
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE INTERACTIONS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Interactions
!
    character*30 :: file_name
    integer      :: io_error
    integer      :: component, component2
!
    file_name='input_interactions'
    open(unit=20,file=file_name,status='old',action='read',iostat=io_error) 
!
    chi = 0.0d0
    if(io_error == 0) then
      read(20,*) kappa
      read(20,*)
      do component = 0, components
        do component2 = component+1, components
          read(20,*) chi(component,component2)
          chi(component2,component) = chi(component,component2)
        end do
      end do
      read(20,*)
      do component = 0, components
        read(20,*) monomer_volume(component)
      end do
    else
      print*, 'Error', io_error, ' while trying to open', file_name 
    end if
!
    do component = 1, components
      do component2 = 1,components
        chi_eff(component2, component) = chi(component2, component) &
           - chi(0,component) - chi(0,component2)
      end do
    end do
    zeta = 1.0d0/monomer_volume(0)
!
  end subroutine Initialize_Interactions
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE INTERACTION ENERGY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  double precision function Interaction_Energy(density)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
!
    double precision, dimension(gridx,gridy,gridz,components) :: volume_fraction
    double precision, dimension(gridx,gridy,gridz)            :: solvent_volume_fraction
    double precision, dimension(gridx,gridy,gridz)            :: local_interaction_energy

    double precision :: field_threshold, energy_threshold
    integer          :: component, component2
!
    field_threshold  = - log(small_number)
    energy_threshold = (small_number*log(small_number)- small_number)
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do

    solvent_volume_fraction = 1.0d0 - sum(volume_fraction,dim=4)
    where (solvent_volume_fraction .ge. small_number)
      local_interaction_energy =  &
        ( solvent_volume_fraction*log(solvent_volume_fraction)- solvent_volume_fraction ) 
    else where
      local_interaction_energy = energy_threshold  &
          - field_threshold * (solvent_volume_fraction - small_number) &
          + 0.5d0*kappa*(solvent_volume_fraction - small_number)**2
    endwhere
    local_interaction_energy = local_interaction_energy * zeta

    do component = 1,components
      local_interaction_energy = local_interaction_energy &
           + chi(0,component) * volume_fraction(:,:,:,component)
      do component2 = 1, components
        local_interaction_energy = local_interaction_energy  &
           + 0.5d0*chi_eff(component,component2) &
              * volume_fraction(:,:,:,component)*volume_fraction(:,:,:,component2)
      end do
    end do

    interaction_energy = dvol*sum(local_interaction_energy)
 
  end function Interaction_Energy
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE MEAN FIELDS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Mean_Fields(density,field)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
    double precision, dimension(gridx,gridy,gridz,components), intent(out):: field
!
    double precision, dimension(gridx,gridy,gridz,components) :: volume_fraction
    double precision, dimension(gridx,gridy,gridz)            :: solvent_volume_fraction
    double precision, dimension(gridx,gridy,gridz)            :: field_0

    double precision :: field_threshold
    integer          :: component, component2
!
    field_threshold = - log(small_number)
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do
!
    solvent_volume_fraction = 1.0d0 - sum(volume_fraction,dim=4)
    where (solvent_volume_fraction .ge. small_number)
      field_0 = - log(solvent_volume_fraction)
    else where
      field_0 = field_threshold - kappa * (solvent_volume_fraction-small_number)
    endwhere
    field_0 = field_0 * zeta
!
    do component = 1,components
      field(:,:,:,component) = field_0 + chi(0,component)
      do component2 = 1, components
        field(:,:,:,component) = field(:,:,:,component) &
          + chi_eff(component,component2)*volume_fraction(:,:,:,component2)
      end do
      field(:,:,:,component) = field(:,:,:,component) &
           * monomer_volume(component)
    end do
 
  end subroutine Calculate_Mean_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Interactions_FloryHuggins_Solution
