!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR FLORY HUGGINS TYPE INTERACTIONS
!               IN A MULTICOMPONENT SYSTEM
!               WITH HELFAND-TYPE INCOMPRESSIBILITY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Interactions_FloryHuggins
!
  use Global
!
  implicit none
!
  private
!
  integer, parameter     :: components = monomer_types
!
  double precision, save :: kappa           ! Helfand compressibility parameter
  double precision, save :: chi(components,components) ! interaction parameters
  double precision, save :: monomer_volume(components) ! monomer volumes
!
  public :: Initialize_Interactions
  public :: Calculate_Mean_Fields
  public :: Calculate_Mean_Ext_Fields
  public :: Interaction_Energy
  public :: Check_Volume_Fraction
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
    if(io_error == 0) then
      read(20,*) kappa
      read(20,*)
      do component = 1, components
        chi(component,component) = 0.0d0
        do component2 = component+1, components
          read(20,*) chi(component,component2)
          chi(component2,component) = chi(component,component2)
        end do
      end do
      read(20,*)
      do component = 1, components
        read(20,*) monomer_volume(component)
      end do
    else
      print*, 'Error', io_error, ' while trying to open', file_name 
    end if
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
    double precision, dimension(gridx,gridy,gridz)            :: total_volume_fraction
    integer   :: component, component2
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do

    total_volume_fraction = sum(volume_fraction,dim=4)
    interaction_energy = 0.5d0 * kappa * dvol*sum((total_volume_fraction-1.0d0)**2)

    do component = 1,components
      do component2 = component+1, components
        interaction_energy = interaction_energy + chi(component,component2) &
          * dvol*sum(volume_fraction(:,:,:,component)*volume_fraction(:,:,:,component2))
      end do
    end do
     print*, 'IE=', interaction_energy
 
  end function Interaction_Energy
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE CONJUGATE FIELDS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Mean_Fields(density,field)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
    double precision, dimension(gridx,gridy,gridz,components), intent(out):: field
!
    double precision, dimension(gridx,gridy,gridz,components) :: volume_fraction
    double precision, dimension(gridx,gridy,gridz)            :: total_volume_fraction
    integer   :: component, component2
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do
!
    total_volume_fraction = sum(volume_fraction,dim=4)
!
    do component = 1,components
!
      field(:,:,:,component) = kappa * (total_volume_fraction-1)
      do component2 = 1, components
        field(:,:,:,component) = field(:,:,:,component) &
          + chi(component,component2)*volume_fraction(:,:,:,component2)
      end do
      field(:,:,:,component) = field(:,:,:,component) * monomer_volume(component)
    end do
 
  end subroutine Calculate_Mean_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CHECK THAT SUM OF AVERAGE VOLUME FRACTIONS GIVES ONE.
!		(FOR USE IN THE CANONICAL ENSEMBLE)
!           THIN FILMS: COUNT ONLY ACCESSIBLE GRID POINTS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Check_Volume_Fraction(density)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
!
    double precision, dimension(gridx,gridy,gridz,components) :: volume_fraction
    double precision  :: average_volume_fraction
    integer           :: component
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do
!
    average_volume_fraction = sum(volume_fraction) / dble(gridx*gridy*gridz)
!
    write(30,*) 'Average total volume fraction is', average_volume_fraction
    print*, 'Average total volume fraction is', average_volume_fraction
    if ( abs(average_volume_fraction - 1.0d0) .gt. 1.0d-3) then
      write(30,*) 'Caution: Differs noticeably from 1'
      print*, 'Caution: Differs noticeably from 1'
    end if
!
  end subroutine Check_Volume_Fraction
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          CALCULATE CONJUGATE FIELD WITH WALL-INTERACTION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine Calculate_Mean_Ext_Fields(density,field)
!
    double precision, dimension(gridx,gridy,gridz,components), intent(in) :: density
    double precision, dimension(gridx,gridy,gridz,components), intent(out):: field
!
    double precision, dimension(gridx,gridy,gridz,components) :: volume_fraction
    double precision, dimension(gridx,gridy,gridz)            :: total_volume_fraction
    integer   :: component, component2
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component) 
    end do
!
    total_volume_fraction = sum(volume_fraction,dim=4)
!
    do component = 1,components
!
      field(:,:,:,component) = kappa*(total_volume_fraction-1)
      do component2 = 1, components
        field(:,:,:,component) = field(:,:,:,component) &
          + chi(component,component2)*volume_fraction(:,:,:,component2)
      end do
      field(:,:,:,component) = field(:,:,:,component) * monomer_volume(component)
    end do

  end subroutine Calculate_Mean_Ext_Fields
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Interactions_FloryHuggins
