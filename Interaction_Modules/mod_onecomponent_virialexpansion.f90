!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR HOMOPOLYMER INTERACTIONS IN SOLUTION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Interactions_OneComponent_VirialExpansion
!
  use Global
!
  implicit none
!
  private
!
  integer, parameter     :: components = 1
  double precision, save :: B    ! Basically the scaled excluded volume
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
!
    if (components.ne.monomer_types) then
      print*, 'Number of monomer types ', monomer_types, &
          & 'does not match number of components in interaction module!'
      print*, 'Interaction module: Interaction_Homopolymers'
      print*, '  Number of components: ', components
      stop
    end if
!    
    file_name='input_interactions'
    open(unit=20,file=file_name,status='old',action='read',iostat=io_error) 
!
    if(io_error == 0) then
      read(20,*) B
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
    interaction_energy = 0.5d0*B*dvol*sum(density*density) 
 
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
    field = B*density
 
  end subroutine Calculate_Mean_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Interactions_OneComponent_VirialExpansion
