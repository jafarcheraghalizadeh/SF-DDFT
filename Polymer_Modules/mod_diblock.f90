!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR DIBLOCK COPOLYMERS 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Molecule_Diblock_Copolymer
!
  use Global
  use Propagate
!
  implicit none
!
  private
  
  integer, parameter  :: block_number = 2
!
  type diblock_copolymer
    integer                                   :: label
    integer, dimension(block_number)          :: monomer_type
    double precision, dimension(block_number) :: block_length
    integer           :: ensemble         ! canonical:1, grand canonical:2
    double precision  :: rho              ! total polymer concentration 
    double precision  :: mu               ! chemical potential of polymer
    double precision  :: energy           ! free energy in SCF field
    integer, dimension(block_number)          :: block_segments
    integer                                   :: segments
    double precision                          :: Qsum
    double precision, dimension(gridx,gridy,gridz,block_number) :: distribution
    double precision, dimension(gridx,gridy,gridz,block_number) :: density
  end type diblock_copolymer
!
  public :: diblock_copolymer
  public :: Initialize_Diblock
  public :: Propagate_Diblock
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZES COPOLYMER
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Diblock(copolymer,label,monomer_type,block_length, &
                                ensemble,rho,mu)
!
   integer, intent(in)                                   :: label, ensemble
   double precision, intent(in)                          :: rho, mu
   double precision, dimension(block_number), intent(in) :: block_length
   integer, dimension(block_number), intent(in)          :: monomer_type
!
   type(diblock_copolymer), intent(out):: copolymer
!
   character*15   :: string_ensemble
   integer :: iAllocStatus
   integer :: block
!
   if (ensemble.eq.1) then
      string_ensemble = 'canonical'
   else if (ensemble.eq.2) then
      string_ensemble = 'grand canonical'
   else
      print*, 'Invalid Ensemble'
   end if
!
   copolymer%label            = label
   copolymer%block_length     = block_length
   copolymer%monomer_type     = monomer_type
   copolymer%ensemble         = ensemble
   copolymer%rho              = rho
   copolymer%mu               = mu
   copolymer%block_segments   = nint(copolymer%block_length/ds)
   copolymer%segments         = sum(copolymer%block_segments)
!
   print*, 'Diblock Copolymer ', copolymer%label, ' is initialized (', &
        &string_ensemble, ' ensemble).'
!
! ----- Initializations for Testing Purposes ---
!
   do block = 1, block_number
     copolymer%density(:,:,:,block) = copolymer%rho &
        * copolymer%block_segments(block)*ds &
        * volume / dvol / dble(gridx*gridy*gridz) ! accounts for boundary effects
   end do
!
  end subroutine Initialize_Diblock
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE PROPAGATORS FOR A DIBLOCK CHAIN OF LENGTH segments*ds
!		DEDUCE SINGLE CHAIN PARTITION FUNCTION Qsum
!               DEDUCE UNNORMALIZED MONOMER DISTRIBUTION
!               DEDUCE MONOMER DENSITY DUE TO THIS CHAIN
!               DEDUCE CHAIN ENERGY OF THIS CHAIN
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Propagate_Diblock(copolymer,expfield)
!
    implicit none

    type(diblock_copolymer), intent(inout)  :: copolymer
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: expfield

    double precision, dimension(:,:,:,:), allocatable   :: qq, qqdag
    double precision, dimension(gridx,gridy,gridz)      :: prop
    integer :: iAllocStatus

    integer :: s,s0,mono,block,monomer_type
!
   allocate(qq(gridx,gridy,gridz,0:copolymer%segments),stat=iAllocStatus)
   allocate(qqdag(gridx,gridy,gridz,0:copolymer%segments),stat=iAllocStatus)

    s0 = 0
    prop=1.0d0
    qq(:,:,:,s0)=prop
    do block = 1,block_number
      monomer_type = copolymer%monomer_type(block)
      do s = s0+1, s0+copolymer%block_segments(block)
         call Propagate_Step(prop,expfield (:,:,:,monomer_type))
         qq(:,:,:,s)=prop
      end do
      s0 = s0 + copolymer%block_segments(block)
    end do
!
    s0 = 0
    prop=1.0d0
    qqdag(:,:,:,s0)=prop
    do block = block_number,1,-1
      monomer_type = copolymer%monomer_type(block)
      do s=s0+1,s0+copolymer%block_segments(block)
         call Propagate_Step(prop,expfield (:,:,:,monomer_type))
         qqdag(:,:,:,s)=prop
      end do
      s0 = s0 + copolymer%block_segments(block)
    end do
!
    mono = copolymer%segments
    copolymer%Qsum = dvol * sum(qq(:,:,:,mono))
!
    s0 = 0
    do block = 1,block_number
      copolymer%distribution(:,:,:,block) = &
          0.5d0*ds*qq(:,:,:,s0)*qqdag(:,:,:,mono-s0)
      do s = s0+1, s0+copolymer%block_segments(block)-1
         copolymer%distribution(:,:,:,block) = copolymer%distribution(:,:,:,block) &
           + ds*qq(:,:,:,s)*qqdag(:,:,:,mono-s)
      end do
      s0 = s0 + copolymer%block_segments(block)
      copolymer%distribution(:,:,:,block) = copolymer%distribution(:,:,:,block) &
         + 0.5d0*ds*qq(:,:,:,s0)*qqdag(:,:,:,mono-s0)
    end do
!
    if (copolymer%ensemble.eq.1) then
      copolymer%density = copolymer%rho*volume * copolymer%distribution/copolymer%Qsum
      copolymer%mu = log(copolymer%rho*volume/copolymer%Qsum)
      copolymer%energy = copolymer%rho*volume * (copolymer%mu - 1.0d0)
    else if (copolymer%ensemble.eq.2) then
      copolymer%density = exp(copolymer%mu) * copolymer%distribution
      copolymer%rho = sum(copolymer%density)*dvol / volume / (ds*copolymer%segments)
      copolymer%energy = - copolymer%rho*volume
    else
      print*, 'Invalid ensemble; Cannot calculate Density'
      stop
    end if
!
    deallocate(qq)
    deallocate(qqdag)
!
    return

  end subroutine Propagate_Diblock
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Molecule_Diblock_Copolymer
