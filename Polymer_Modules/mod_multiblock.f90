!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR LINEAR MULTIBLOCK COPOLYMERS 
!
!           6.7.2020
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Molecule_Multiblock_Copolymer
!
  use Global
  use Propagate
!
  implicit none
!
  private
!
  type multiblock_copolymer
    integer                                    :: label
    integer                                    :: block_number
    integer, dimension(:),allocatable          :: monomer_type
    double precision, dimension(:),allocatable :: block_length
    integer           :: ensemble         ! canonical:1, grand canonical:2
    double precision  :: rho              ! total polymer concentration 
    double precision  :: mu               ! chemical potential of polymer
    double precision  :: energy           ! free energy in SCF field
    integer, dimension(:),allocatable   :: block_segments
    integer           :: segments
    double precision  :: Qsum
    double precision, dimension(:,:,:,:),allocatable :: distribution
    double precision, dimension(:,:,:,:),allocatable :: density
  end type multiblock_copolymer
!
  public :: multiblock_copolymer
  public :: Initialize_Multiblock
  public :: Propagate_Multiblock
!
contains

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZES COPOLYMER
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Multiblock(copolymer,label,nblocks, &
                    monomer_type,block_length, ensemble,rho,mu)
!
   implicit none

   integer, intent(in)                 :: label, nblocks, ensemble
   double precision, intent(in)        :: rho, mu
   double precision, dimension(:), intent(in) :: block_length
   integer, dimension(:), intent(in)          :: monomer_type
!
   type(multiblock_copolymer), intent(out):: copolymer
!
   character*15   :: string_ensemble
   integer :: iAllocStatus
   integer :: block
!
   allocate(copolymer%monomer_type(nblocks),stat=iAllocStatus)
   allocate(copolymer%block_length(nblocks),stat=iAllocStatus)
   allocate(copolymer%block_segments(nblocks),stat=iAllocStatus)
   allocate(copolymer%distribution(gridx,gridy,gridz,nblocks),stat=iAllocStatus)
   allocate(copolymer%density(gridx,gridy,gridz,nblocks),stat=iAllocStatus)
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
   copolymer%block_number     = nblocks
   copolymer%monomer_type     = monomer_type
   copolymer%block_length     = block_length
   copolymer%ensemble         = ensemble
   copolymer%rho              = rho
   copolymer%mu               = mu
   copolymer%block_segments   = nint(copolymer%block_length/ds)
   copolymer%segments         = sum(copolymer%block_segments)
!
   print*, 'Multiblock Copolymer ', copolymer%label, ' is initialized (', &
        &string_ensemble, ' ensemble).'
!
! ----- Initializations for Testing Purposes ---
!
   do block = 1, nblocks
     copolymer%density(:,:,:,block) = copolymer%rho &
        * copolymer%block_segments(block)*ds &
        * volume / dvol / dble(gridx*gridy*gridz) ! accounts for boundary effects
   end do
!
  end subroutine Initialize_Multiblock
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE PROPAGATORS FOR A DIBLOCK CHAIN OF LENGTH segments*ds
!		        DEDUCE SINGLE CHAIN PARTITION FUNCTION Qsum
!               DEDUCE UNNORMALIZED MONOMER DISTRIBUTION
!               DEDUCE MONOMER DENSITY DUE TO THIS CHAIN
!               DEDUCE CHAIN ENERGY OF THIS CHAIN
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
  subroutine Propagate_Multiblock(copolymer,expfield)
     use, intrinsic :: ieee_arithmetic
    implicit none

    type(multiblock_copolymer), intent(inout)  :: copolymer
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: expfield

    double precision, dimension(:,:,:,:), allocatable   :: qq, qqdag
    double precision, dimension(gridx,gridy,gridz)      :: prop
    integer :: iAllocStatus, nblocks
    double precision::hu
    integer :: s,s0,mono,block,monomer_type
logical:: inf
 double precision::ma=0.0
    allocate(qq(gridx,gridy,gridz,0:copolymer%segments),stat=iAllocStatus)
    allocate(qqdag(gridx,gridy,gridz,0:copolymer%segments),stat=iAllocStatus)
!
    nblocks = copolymer%block_number
!
inf =.True.
    s0 = 0
    prop=1.0d0
    qq(:,:,:,s0)=prop
    do block = 1,nblocks
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
    do block = nblocks,1,-1
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
   if (copolymer%Qsum>ma)then
   ma=copolymer%Qsum
   end if 
!
    s0 = 0
    do block = 1,nblocks
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

    deallocate(qq)
    deallocate(qqdag)
!
    if (copolymer%ensemble.eq.1) then
      copolymer%density = (copolymer%rho*volume* copolymer%distribution/copolymer%Qsum)
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
    return

  end subroutine Propagate_Multiblock
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Molecule_Multiblock_Copolymer
