!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM (AB HOMOPOLYMER BLEND)
!          
!           DYNAMIC DENSITY FUNCTIONAL CALCULATION
!           Reference (chain dynamics): Qi, Schmid, Macromolecules 50, 9831 (2017)
!                https://doi.org/10.1021/acs.macromol.7b02017
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_blend
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use Global
  use Molecule_Homopolymer
  use Interactions_FloryHuggins
  use mtmod
!
  implicit none
!
! GLOBAL PARAMETERS: DECLARED / DEFINED IN MODULE GLOBAL
!
! double precision, parameter :: pi = 3.141592653589793d0
!
! integer, parameter          :: monomer_types=2
!
! integer, parameter  :: gridx, gridy, gridz ! lattice points in each direction
! double precision    :: sizex, sizey, sizez ! system size
! double precision    :: volume	             ! system volume
! double precision    :: dvol 		     ! volume element
! double precision    :: ds		     ! discretization along chain
!
  integer, parameter  :: polymer_types=monomer_types

  type(homopolymer), dimension(polymer_types) :: polymer
 
  double precision :: Ctotal    ! total dimensionless concentration of polymers

  double precision, dimension(polymer_types) :: mu  ! chemical potential of polymer 
  double precision, dimension(polymer_types) :: C   ! dimensionless concentration of polymer 
  double precision, dimension(polymer_types) :: number_of_chains
  double precision, dimension(polymer_types) :: chain_length

  ! double precision    ::  energy, interaction_energy
  double precision    ::  energy

  integer             ::  number_of_timesteps ! number of timesteps
  double precision    ::  timestep            ! length of timestep
  double precision    ::  wcutoff        ! cutoff for DDFT calculations (numerical stability)
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: density
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: field
  double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: mobility_qq
  double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: mobility
  double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types),save :: lambda_a
  logical          :: qexist, qread, qepd,qdmin
  character*24     :: file_name
  integer          :: polymer_label, monomer_type,gcfstat
  integer          :: polymer_ensemble    ! canonical:1, grand canonical:2
  integer          :: qinit   ! 1:sharp interface, 2: random densities, 3: random fields
  double precision :: dx, dy, dz
  double precision, dimension(gridx/2+1,gridy,gridz) :: q1
  integer          :: io_error, x,y,z,step, rgx,rgy,rgz, i
  real::r
  character(len=20) :: filenameq,filenameq_dens

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     INITIALIZATIONS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! ------ READ INPUT PARAMETERS 
 call get_command_argument(1, filenameq)
 call get_command_argument(2, filenameq_dens)
open(unit=200,file=trim(filenameq),action='write',iostat=io_error)
open(unit=100,file='E_'//trim(filenameq),action='write',iostat=io_error)
 call init_random_seed()
  file_name='input_general'
  call read_input_parameters(file_name)
!
! ------ GENERAL INITIALIZATIONS 
!
  call Initialize

  C = number_of_chains/volume
  Ctotal = sum(C)

  dx = sizex/dble(gridx)
  dy = sizey/dble(gridy)
  dz = sizez/dble(gridz)
!
! ------ INITIALIZE HOMOPOLYMERS
!
  do polymer_label = 1,2
    monomer_type = polymer_label
    polymer_ensemble  = 1  ! canonical
!   polymer_ensemble  = 2  ! grand canonical
    call Initialize_Homopolymer(polymer(monomer_type), polymer_label, &
                     monomer_type, chain_length(polymer_label), &
                     polymer_ensemble, C(polymer_label), mu(polymer_label))
  end do
!
! ------ INITIALIZE INTERACTIONS
!
  call Initialize_Interactions
!
 if (qepd) then
    !file_name='input_mobilities_d'
   file_name='Mobilityentangledprint'
    call Read_Mobilities(file_name)
    call calculate_lambda(lambda_a,mobility)
  end if
! ------ PREPARE OUTPUT FILES
!
  inquire(file='check.dat', exist=qexist)
  if (qexist) then
    open(unit=30,file='check.dat', status='old',action='write',position='append',iostat=io_error)
  else
    open(unit=30,file='check.dat', status='new',action='write',iostat=io_error)
  end if
  if(io_error /=0) then                   
        write(*,*) 'Error', io_error , ' while trying to open check.dat'
  end if

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   INITIAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  rgx = nint(chain_length(1)/dx)
  rgy = nint(chain_length(1)/dy)
  rgz = nint(chain_length(1)/dz)
  if (qread) then
    open (unit=33, file='cnf.in', status='old', action='read', iostat=io_error)
    if(io_error /=0) then                   
      write(*,*) 'Error', io_error , ' while trying to open cnf.in'
    end if
    read (33,*) field
    close (unit=33)
  else
    if (qinit.eq.1) then
       do z=1,gridz
         do y=1,gridy
           do x=1,gridx/2+1
             density(x,y,z,1) = 1.
             density(x,y,z,2) = 0.
           end do
           do x=gridx/2+1,gridx
             density(x,y,z,1) = 0.
             density(x,y,z,2) = 1.
           end do
         end do
       end do
       call Get_Conjugate_Fields(density,field,0.01d0,10000,3,gcfstat)
     else if (qinit .eq. 2) then
       do i = 1,monomer_types
         do z=1,gridz
           do y=1,gridy
              do x=1,gridx
                density(x,y,z,i) = grnd()     
              end do
           end do
         end do
       end do
      density = Ctotal*(1.0d0 + (2.0d0*density -1.0d0)*0.001d0)
      call Get_Conjugate_Fields(density,field,0.01d0,10000,3,gcfstat)
     elseif (qinit .eq. 3) then
       do i = 1,monomer_types
         do z=1,gridz
           do y=1,gridy
              do x=1,gridx
              call random_number(r)
                field(x,y,z,i) = r!grnd()     
              end do
           end do
         end do
       end do

 density=0.0d0
     call Calculate_Densities(field,density)
     else if (qinit .eq. 4) then
     open (unit=165, file=trim(filenameq_dens), status='old', action='read', iostat=io_error)
     if(io_error /=0) then
          write(*,*) 'Error', io_error , ' while trying to open dens.in'
          print*,filenameq_dens
          end if
          do x=1,gridx
            do y=1,gridy
               do z=1,gridz
                 read (165,*) density(x,y,z,1), density(x,y,z,2)
                 !print*,density(x,y,z,1),density(x,y,z,2)
                end do
            end do
          end do
        print*, 'Successful reading of density'
        close (unit=165)
        call Get_Conjugate_Fields(density,field,0.01d0,2000,3,gcfstat)
        call Calculate_Densities(field,density)
        call Get_Conjugate_Fields(density,field,0.01d0,2000,3,gcfstat)

     else
       print*, 'Initialization not well-defined!'
       stop
     end if
   end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   EXTERNAL POTENTIAL DYNAMICS SIMULATION / SCF CALCULATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  if (qepd) then

    write(30,*) 'BEGIN DYNAMIC SCF SIMULATION'
    write(30,*) 
    qdmin = .true.
    do step = 1,number_of_timesteps
       call DDFT_Timestep_Euler(density,field,lambda_a,timestep,qdmin)
        
     !  call Chain_Dynamics_Timestep_Euler(field,density,timestep)
     ! call EPD_Timestep_Euler(field,timestep)
     ! call EPD_Timestep_RungeKutta(field,timestep)
     ! call Calculate_Densities(field,density)
      if (((step/50)*50).eq. step) then                  ! FOR TESTING
       call Calculate_Energy                             ! FOR TESTING
       print*, step, energy ,lambda    ! FOR TESTING
  do z=1,gridz
           do y=1,gridy
              do x=1,gridx
                 write(200, *) density(x,y,z,1),density(x,y,z,2)!,field(x,y,z,1),field(x,y,z,2)
              end do

           end do
        end do
         write(100,*)timestep*step,energy
     end if                                             ! FOR TESTING

    end do

 else

    write(30,*) 'BEGIN SCF CALCULATION'
    write(30,*) 

    call Find_SCF_Saddle_Point(density,field)
    call Calculate_Densities(field,density)

  end if

  call Calculate_Energy

  if (qepd) then
    write(30,*) '   Final values after ', number_of_timesteps, 'steps'
    print*, ' Final values after ', number_of_timesteps, 'steps'
  else
    write(30,*) '   Final values '
    print*, ' Final values '

  end if
  write(30,*) '   Energy: ', energy, ' C: ', polymer%rho, ' mu_scf: ', polymer%mu
  write(30,*) 
!
  print*, '   Energy: ', energy, ' C: ', polymer%rho, ' mu_scf: ', polymer%mu
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF FINAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  call Calculate_Densities(field,density)

      open(unit=60,file='dens.dat', status='replace',action='write',iostat=io_error)
      open(unit=61,file='field.dat', status='replace',action='write',iostat=io_error)
        do z=1,gridz
           do y=1,gridy
              do x=1,gridx
                 write(60, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', (density(x,y,z,i),i=1,monomer_types)
                 write(61, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', (field(x,y,z,i),i=1,monomer_types)
              end do
             write(60, *) 
             write(61, *) 
           end do
        end do
      close(unit=60)
      close(unit=61)
!
      open(unit=33, file='cnf.out', status = 'replace', action='write', iostat=io_error)
        write (33,*) field
      close (unit=33)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   FINISH (CLOSE FILES ETC.)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  close (unit=30)
  call finish
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   END OF MAIN PROGRAM. NOW INNER SUBROUTINES
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO READ INPUT PARAMETERS  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine read_input_parameters(file_name)
    
    character*24, intent(in) :: file_name
    integer                  :: io_error

    open(unit=20,file=file_name,status='old',action='read',iostat=io_error)

      if(io_error ==0) then

        read(20,*) sizex
        read(20,*) sizey
        read(20,*) sizez
        read(20,*)
        do polymer_label = 1,polymer_types
          read(20,*) number_of_chains(polymer_label)
          read(20,*) mu (polymer_label)
          read(20,*) chain_length(polymer_label)
          read(20,*)
        end do
        read(20,*) ds
        read(20,*)
        read(20,*) number_of_timesteps
        read(20,*) timestep
        read(20,*)
        read(20,*) qread  ! qread = .true. : read from file cnf.in
        read(20,*) qinit  ! only for qread=.false.  
            ! 1=interface, 2=random densities, 3=random fields
        read(20,*) qepd   ! qepd = .true. : EPD simulation; otherweise SCF calculation
        read(20,*) wcutoff
      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if

    close(unit=20)
  
    return

  end subroutine read_input_parameters
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO CALCULATE DENSITIES  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Densities(field,density)
!
   use Propagate
implicit none
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: expfield

    integer :: monomer_type, i

    call Calculate_Expfield(field,expfield)
    do i = 1,polymer_types
      call Propagate_Homopolymer(polymer(i),expfield)
    end do
 
    density = 0.0d0
    do i = 1,polymer_types
       monomer_type = polymer(i)%monomer_type
       density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + polymer(i)%density(:,:,:)
    end do
 
    return

  end subroutine Calculate_Densities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO CALCULATE ENERGY  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Energy
!
    double precision :: chain_energy
    integer :: i

    chain_energy = - dvol*sum(density*field)
    do i = 1,polymer_types
       chain_energy = chain_energy + polymer(i)%energy 
    end do
    energy = chain_energy + interaction_energy(density)                

    return

  end subroutine Calculate_Energy
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO CALCULATE FIELDS CORRESPONDING TO GIVEN DENSITIES  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Get_Conjugate_Fields(density,field,accuracy_scf,itermax_scf,mixing_type,gcfstat)
!
   use Propagate
   use global
   use iterate
     implicit none
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)    :: density
    integer, intent(out) :: gcfstat
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity

    double precision ::  accuracy, accuracy_scf
    integer :: mixing_type, iter_status, iter, itermax_scf

    !accuracy_scf = 0.01d0
    !mixing_type = 3
    !iter_max = 1000
    !lambda = 0.5

!    small = 1.d-6
!    where (density > small)
!      density0 = density
!    elsewhere
!      density0 = small
!    end where

!lambda=0.5
  iter_status = 1
    do iter = 1,itermax_scf

 if (norm(2)-norm(1)>0 .and. iter>100)then
      lambda=max(lambda_min,0.9d0*lambda)
      !print*,'i am here'
  else
      lambda=min(lambda_max,1.1d0*lambda)
     ! print*,lambda_max
  end if



      call Calculate_Densities(field,newdensity)
      dfield = newdensity - density

      call Check_Accuracy(dfield,accuracy)
!      if ( ( (iter/100)*100) .eq. iter) then
!        print*, iter, accuracy
!      end if

      if (.not.(accuracy > accuracy_scf*0.5) ) then  ! double negation to deal with NaNs
!        PRINT*, 'exit Get_Conjugate_Field after', iter, 'steps '
!        PRINT*, 'Accuracy: ', accuracy
         exit
       end if
      
       iter_status = min(iter,mixing_dim+1)
       call Mix_Fields(mixing_type,field,dfield,iter_status)

     end do

     if (.not.(accuracy < accuracy_scf*0.5)) then  ! double negation to deal with NaNs
       if (.not.(accuracy <  accuracy_scf)) then  ! double negation to deal with NaNs
       !  print*, 'Warning: Conjugate fields not found within ',iter,' iteration steps!'
      !   print*, 'Accuracy: ', accuracy
         gcfstat = 2 ! iteration not successful
         if (.not.(accuracy <  200.0d0)) then  ! double negation to deal with NaNs
           stop
         end if
       else
       !  print*, 'Caution: Finding conjugate fields takes long!'
        ! print*, 'Accuracy after ', iter, 'iteration steps: ', accuracy
         gcfstat = 1 ! iteration slow
       end if
     else
       gcfstat = 0 ! iteration successful
!       print*, 'Conjugate fields found within ',iter,' iteration steps!'
!       print*, 'Accuracy: ', accuracy
     end if

    return

  end subroutine Get_Conjugate_Fields
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO FIND SCF SADDLE POINT
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Find_SCF_Saddle_Point(density,field)
!
   use Propagate
!   use Interactions_FloryHuggins

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity, newfield

    double precision :: accuracy, accuracy_scf
    integer :: mixing_type, iter_status, iter, iter_max

    accuracy_scf = 0.0000001
    mixing_type = 2
    iter_max = 2000
    !lambda = 0.2

    iter_status = 1

    do iter = 1,iter_max

      call Calculate_Densities(field,newdensity)
      call Calculate_Mean_Fields(newdensity,newfield)
      dfield = newfield - field
      call Check_Accuracy(dfield,accuracy)
!      if ( ( (iter/100)*100) .eq. iter) then
!        print*, iter, accuracy
!      end if

      if (.not.(accuracy > accuracy_scf) ) then  ! double negation to deal with NaNs
        print*, 'SCF saddle point found after ',iter,' iteration steps!'
        print*, 'Accuracy: ', accuracy
        exit
      end if
      call Mix_Fields(mixing_type,field,dfield,iter_status)
 
    end do
    
    if (.not.(accuracy < accuracy_scf)) then  ! double negation to deal with NaNs
      print*, 'Caution: SCF saddle point not found within ',iter,' iteration steps!'
      print*, 'Accuracy: ', accuracy
    end if

    density = newdensity
 
    return

  end subroutine Find_SCF_Saddle_Point
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO MIX FIELDS IN SCF ITERATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE ACHIEVED ACCURACY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Check_Accuracy(dfield,accuracy)

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: dfield
    double precision, intent(out)       :: accuracy

    accuracy=sum(abs(dfield))*dvol

    return
  end subroutine check_accuracy
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   EPD TIME STEP 
!          USING FORWARD EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   MODIFIED EPD TIME STEP 
!          USING FORWARD EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  
  
subroutine Read_Mobilities(file_name)
    
    use Fourier_fftw3

    implicit none

    character*24, intent(in) :: file_name

    integer, parameter       :: nqmax=1000
    double precision, dimension(nqmax,monomer_types,monomer_types) :: mobility_table
    double precision,allocatable::q_values(:)
    double precision :: dum, dq, qq, q, small
    integer          :: io_error, nq, iq, x, y, z,type1,type2,idx
    double precision::q_low,q_high,d_low,d_high
!   ------------ read mobilities 
!
    small = 1.d-6
    open(unit=150,file='mob.dat')
    open(unit=20,file=file_name,status='old',action='read',iostat=io_error)
      if(io_error ==0) then
        read(20,*) nq
        allocate(q_values(nq))
          if (nq.gt.nqmax) then
            stop 'too many entries in mobility file -> adjust nqmax'
          end if
        read(20,*) dq
        read(20,*)
        read(20,*)
        read(20,*)
        do iq = 1, nq
           read(20,*) dum , mobility_table(iq,1,1)  &
                     ,mobility_table(iq,1,2), mobility_table(iq,2,2)
          q_values(iq)=dum
           if (abs(dum -iq*dq) .ge. small ) then
             print*, 'q values in mobility are not consistent', dum, iq*dq,dq
             stop 
           end if
        end do
        mobility_table(:,2,1) = mobility_table(:,1,2)
      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if
    close(unit=20)
!
!   --------------- Calculate mobility_qq = Lambda(q)*q^2
!
  ! do type1=1,monomer_types
     !    do type2=1,monomer_types
    !do z = 1,gridz
    ! do y = 1,gridy
     ! do x = 1,gridx/2+1
        ! qq = real(-laplace(x,y,z), kind(0.0d0))
        ! q  = dsqrt(qq)
       !  iq = int(q/dq)
        ! if (iq .eq. 0) then
        !   mobility_qq(x,y,z,type1,type2) = mobility_table(1,type1,type2)*q/dq
        ! else if (iq .ge. nq) then
           !mobility_qq(x,y,z,type1,type2) = mobility_table(nq,type1,type2)
        ! else 
          ! mobility_qq(x,y,z,type1,type2)  &
             ! = (   mobility_table(iq,type1,type2)  *(q  -  iq*dq) &
               !   + mobility_table(iq+1,type1,type2)*((iq+1)*dq-q) ) / dq
        ! end if
       ! if ((mobility_qq(x,y,z,type1,type2)>=0))then
        !mobility_qq(x,y,z,type1,type2) = min(qq*mobility_qq(x,y,z,type1,type2),wcutoff)
        !else if ((mobility_qq(x,y,z,type1,type2)<0))then
        !mobility_qq(x,y,z,type1,type2) = max(qq*mobility_qq(x,y,z,type1,type2),-wcutoff+1.0e-15)
        ! endif
      !end do
     !end do
    !end do
   ! enddo
   ! enddo
       !=============================================
    do type1=1,monomer_types
        do type2=1,monomer_types
          do z = 1,gridz
            do y = 1,gridy
               do x = 1,gridx/2+1
                  qq = real(-laplace(x,y,z), kind(0.0d0))
                  q1(x,y,z)=dsqrt(qq)
                  q  = dsqrt(qq)
                  iq = int(q/dq)
                  idx = 1
                  do while (idx <= nq)
                      if (q_values(idx) > q) exit
                       idx = idx + 1
                  end do
                    idx = idx - 1
                 if (idx <= 1) then
                    mobility(x,y,z,type2,type1) = mobility_table(1,type1,type2)
                  elseif (idx >= nq) then
                      mobility(x,y,z,type2,type1) = mobility_table(nq,type1,type2)
                  else
                      q_low  = q_values(idx)
                      q_high = q_values(idx + 1)
                      d_low  = mobility_table(idx,type1,type2)
                      d_high = mobility_table(idx + 1,type1,type2)
                      ! Perform linear interpolation
                      mobility(x,y,z,type1,type2) = d_low + (q - q_low) * (d_high - d_low) / (q_high - q_low)
                      if ((q_high - q_low)==0)print*,'ERROR Mobility: Zero in denominator'
                  end if
            
            if ((mobility(x,y,z,type1,type2)>=0))then
        mobility_qq(x,y,z,type1,type2) = min(qq*mobility(x,y,z,type1,type2),wcutoff)
        else if ((mobility(x,y,z,type1,type2)<0))then
        mobility_qq(x,y,z,type1,type2) = max(qq*mobility(x,y,z,type1,type2),-wcutoff+1.0e-15)
         endif
             !    if (min( mobility(x,y,z,type2,type1) * q**2, wcutoff)==0.0)print*,x,y,z

          

     end do
     end do
    end do
    enddo
    enddo
 end subroutine Read_Mobilities

 subroutine calculate_lambda(lambda_a,mobility)
    implicit none
   double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types),intent(in) ::mobility
   double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types),intent(out) ::lambda_a
   !double precision, dimension(gridx,gridy,gridz)   :: dummy
   !double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier
   !  double precision, dimension(:,:,:,:,:),intent(in) ::mobility
    integer :: x, y, z

    ! Initialize lambda_a to zero
    lambda_a(:,:,:,:,:) = 0.0d0

    ! Calculate lambda_a(:,:,:,1,1)
 

    ! Apply conditions where lambda_a(:,:,:,1,1) is not zero
    do x = 1, gridx/2+1
        do y = 1, gridy
            do z = 1, gridz
                 lambda_a(x,y,z,1,1) = dsqrt(mobility(x,y,z,1,1))
                if (lambda_a(x,y,z,1,1) /= 0.0d0) then
                    ! Calculate lambda_a(:,:,:,2,1)
                    lambda_a(x,y,z,2,1) = mobility(x,y,z,2,1) / lambda_a(x,y,z,1,1)
                else
                print*, 'Error: lambda_a(:,:,:,1,1) is zero'
                end if
                    ! Calculate lambda_a(:,:,:,2,2)
                    lambda_a(x,y,z,2,2) = dsqrt(mobility(x,y,z,2,2) - lambda_a(x,y,z,2,1)**2)
            end do
        end do
    end do


end subroutine calculate_lambda
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      DDFT TIME STEP
!          USING FORWARD EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   subroutine DDFT_Timestep_Euler (density,field,lambda_a,timestep,qdmin)
 !
    use Fourier_fftw3

      implicit none
      integer,parameter:: Dimention =3
      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density
      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field

      !double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types), intent(in) :: mobility_qq
      !double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types), intent(in) :: mobility
      double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types), intent(in) :: lambda_a
      double precision, intent(in) :: timestep
      logical, intent(in)          :: qdmin

      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: newfield
      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: mulocal
      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: ddensity
      double precision, dimension(gridx,gridy,gridz,monomer_types,Dimention)   :: h_r
      !double precision, dimension(gridx,gridy,gridz,monomer_types)   :: H_k_gamma


      double complex, dimension(gridx/2+1,gridy,gridz,monomer_types,Dimention) :: mu_q,f_k
      double complex, dimension(gridx/2+1,gridy,gridz,monomer_types,dimention) :: H_k
      double complex, dimension(gridx/2+1,gridy,gridz,monomer_types) :: ddensity_q,H_k_gamma,density_qq,j_a_qq
      double complex, dimension(gridx/2+1,gridy,gridz,monomer_types,dimention) :: j_a
      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: j_a_r
         double precision, dimension(gridx,gridy,gridz,monomer_types,dimention)   :: j_a_r_v
      double precision, dimension(gridx,gridy,gridz)   :: dummy
      double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier,dummy_c

!
      double precision :: dsmall
      integer :: monomer_type,type1,type2, gcfstat
      double complex,parameter::i1=(0.0d0,1.0d0)  !! sqrt(-1)  = i in Complex number
      double complex, dimension(gridx/2+1,gridy, gridz,monomer_types,3) :: rand_F
      integer::dim
!     Local potential in Real space

      call Get_Conjugate_Fields(density,field,0.01d0,2000,3,gcfstat)
      call Calculate_Mean_Fields(density,newfield)
       mulocal = newfield - field
!      do monomer_type = 1,monomer_types
!        offset(monomer_type) = sum(mulocal(:,:,:,monomer_type))/float(ngrid)
!        mulocal(:,:,:,monomer_type) = mulocal(:,:,:,monomer_type) &
!             - offset(monomer_type)
!      end do

!     Local potential in Fourier space
!=====================================================================1-2
do dim=1,Dimention   
    do monomer_type = 1, monomer_types
        dummy = mulocal(:,:,:,monomer_type)
        call Real2Fourier(dummy,dummyfourier)
        dummyfourier = dummyfourier*nabla(:,:,:,dim)  ! g_b (k)= i* k* mu_b (k) (1)
        call Fourier2Real(dummy,dummyfourier)
        dummy = dummy*dsqrt(abs(density(:,:,:,monomer_type)))
        call Real2Fourier(dummy,dummyfourier)
        mu_q(:,:,:,monomer_type,dim)= dummyfourier !g_b (k')
  end do
end do

!=========================================================================2
!do monomer_type = 1, monomer_types
      !  dummy = dsqrt(density(:,:,:,monomer_type))
        !call Real2Fourier(dummy,dummyfourier)
        !call real2fourier_complex(dummyfourier, dummy_c)
       ! density_qq(:,:,:,monomer_type)= dummy_c  !  rho(k')
!enddo

!==========================================================================3
!do dim=1,Dimention
  !do  monomer_type = 1, monomer_types
  !  dummyfourier= mu_q(:,:,:,monomer_type,dim)*density_qq(:,:,:,monomer_type)
!call Fourier2Real_complex(dummy_c,dummyfourier)
!f_k(:,:,:,monomer_type,dim)=dummy_c   !F^{-1}  g(k') *sqrt(rho(k'))
!  end do
!end do
!print*,ddensity_q(5,5,1,1),f_k(5,5,1,1,1)

!===========================================================================4
     h_k(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
do dim=1,Dimention
  do  type1 = 1, monomer_types
      do  type2 = 1, monomer_types
      h_k(:,:,:,type1,dim) = h_k(:,:,:,type1,dim) + mobility(:,:,:,type2,type1)*mu_q(:,:,:,type2,dim)
      end do
  
  dummyfourier=h_k(:,:,:,type1,dim)
  !mu_q1(:,:,:,monomer_type,dim)=dummyfourier
  call Fourier2Real(dummy,dummyfourier)
  h_r(:,:,:,type1,dim)=dummy
  enddo
  enddo

!==========================================================================5
    H_k_gamma(:,:,:,:)=dcmplx(0.0d0,0.0d0)
do  monomer_type = 1, monomer_types
  do dim=1,Dimention
       h_r(:,:,:,monomer_type,dim)=h_r(:,:,:,monomer_type,dim)*dsqrt(abs(density(:,:,:,monomer_type)))
       dummy=h_r(:,:,:,monomer_type,dim)
       call Real2Fourier(dummy,dummyfourier)
       H_k_gamma(:,:,:,monomer_type)=H_k_gamma(:,:,:,monomer_type)+dummyfourier*nabla(:,:,:,dim)
  end do
end do
!print*, H_k(5,5,1,1)
!print*,mu_q1(5,5,1,1,1),h_r(5,5,1,1,1)

!!==========================================================================
     ddensity_q = 0.0d0
     !lambda_a(:,:,:,:,:)=0.0d0
     !lambda_a(:,:,:,1,1)=dsqrt((mobility(:,:,:,1,1)))
     !WHERE (lambda_a(:,:,:,1,1) /= 0.0)
     !lambda_a(:,:,:,2,1) = mobility(:,:,:,2,1) / lambda_a(:,:,:,1,1)
     !END WHERE
    ! print*,sqrt((mobility(:,:,:,2,2)- lambda_a(x,y,z,2,1)**2))
    ! lambda_a(:,:,:,2,2)=dsqrt((mobility(:,:,:,2,2)- lambda_a(:,:,:,2,1)**2))
!!==========================================================================
 	  ! WHERE ((mobility(:,:,:,2,2)- lambda_a(:,:,:,2,1)**2)< 0.0)
   !  stop 1
   !  END WHERE
!print*,ddensity_q(5,5,1,1,1),'--'
   

call random_stdnormal(rand_F,KT,0.0d0)
 	!call random_stdnormal(rand_F1,kt,0.0d0)

 	   !do z = 1,gridz
            !  do y = 1,gridy
 	  		! do x = 1,(gridx/2)+1
        j_a=dcmplx(0.0d0,0.0d0)
do dim=1,dimention
	! J_1=============================================
     !	j_a(x,y,z,1)=j_a(x,y,z,1)&
     !	+(rand_F(x,y,z,1))*nabla(x,y,z,1)* lambda_a(x,y,z,1,1)+(rand_F(x,y,z,2))&
    ! 	*nabla(x,y,z,2)* lambda_a(x,y,z,1,1)!+(rand_F(x,y,z,3))*nabla(x,y,z,3)*lambda_a(x,y,z,1,1)
j_a(:,:,:,1,dim)=j_a(:,:,:,1,dim)+ (lambda_a(:,:,:,1,1))*rand_F(:,:,:,1,dim)!*sqrt(ddensity_q(x,y,z,1,1))
     	
        ! J_2===================gridx=========================
 	!j_a(x,y,z,2)=j_a(x,y,z,2)&
 !	+(rand_F(x,y,z,1))*nabla(x,y,z,1)* lambda_a(x,y,z,2,1)+(rand_F(x,y,z,2))&
 	!*nabla(x,y,z,2)* lambda_a(x,y,z,2,1)!+(rand_F(x,y,z,3))*nabla(x,y,z,3)* lambda_a(x,y,z,2,1))
j_a(:,:,:,2,dim)=j_a(:,:,:,2,dim)+ (lambda_a(:,:,:,2,1))*rand_F(:,:,:,1,dim)
	
     !	j_a(x,y,z,2)=j_a(x,y,z,2)&
     !!	+(rand_F1(x,y,z,1))*nabla(x,y,z,1)* lambda_a(x,y,z,2,2)+(rand_F1(x,y,z,2))&
     !	*nabla(x,y,z,2)* lambda_a(x,y,z,2,2)!+(rand_F1(x,y,z,3))*nabla(x,y,z,3)*lambda_a(x,y,z,2,2))
     		
j_a(:,:,:,2,dim)=j_a(:,:,:,2,dim)+ (lambda_a(:,:,:,2,2))*rand_F(:,:,:,2,dim)
      enddo

j_a_qq=dcmplx(0.0d0,0.0d0)
do monomer_type=1,monomer_types
    do dim = 1,dimention
          dummyfourier=j_a(:,:,:,monomer_type,dim)
          dummy=0.0d0
          call Fourier2Real(dummy,dummyfourier)
          dummy=dummy*dsqrt(abs(density(:,:,:,monomer_type)))
          call Real2Fourier(dummy,dummyfourier)
          j_a_qq(:,:,:,monomer_type)= j_a_qq(:,:,:,monomer_type)+dummyfourier*nabla(:,:,:,dim)
    enddo
enddo

     ! j_a(x,y,z,1)=j_a(x,y,z,1)*i1
    !   j_a(x,y,z,2)=j_a(x,y,z,2)*i1
  ! j_a=j_a*i1
   !print*,j_a
!     Density derivative in Real space
 j_a_r=0.0d0
do monomer_type = 1, monomer_types
    dummyfourier = H_k_gamma(:,:,:,monomer_type)
    call Fourier2Real(dummy,dummyfourier)
    ddensity(:,:,:,monomer_type) = dummy
    dummyfourier = j_a_qq(:,:,:,monomer_type)
    call Fourier2Real(dummy,dummyfourier)
    j_a_r(:,:,:,monomer_type)=sqrt(timestep)*dummy*sqrt(2.0d0)

end do
        
   
        ! print*,  h_k(5,5,1,1)
   	 !  print*, j_a(7,7,1,2)
!     Euler step
dsmall = 1.d-4
     !print*,j_a_r(7,7,1,1)

do monomer_type = 1, monomer_types
    
    !  call random_stdnormal_real(ran,1.0d0,0.0d0,gridx,gridy,gridz)
    density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + &
    timestep*(ddensity(:,:,:,monomer_type))+j_a_r(:,:,:,monomer_type)
     !sqrt(timestep)*ran*0.1
enddo
     !print*,density(5,5,1,1),density(5,5,1,2)
     ! dummy = density(:,:,:,1)
      !  call Real2Fourier(dummy,dummyfourier)
        
           !  do z = 1,gridz
            !  do y = 1,gridy
 	  	!	 do x = 1,gridx
         !write(250,*) real(dummyfourier(x,y,z)), aimag(dummyfourier(x,y,z))
         !  write(250,*)density(x,y,z,1),density(x,y,z,2)
       ! enddo
        ! enddo
        ! enddo
   !   if (qdmin) then
      !  dum=minval(density)
      !  if (dum .le. dsmall) then
         ! density = max(density, dsmall)
       ! end if
    !  end if
call Get_Conjugate_Fields(density,field,0.01d0,1000,3,gcfstat)

if (qdmin .and. (gcfstat.ge.1)) then
  call Calculate_Densities(field,density) 
end if
!
return
 
end subroutine DDFT_Timestep_Euler
   
!
  subroutine Mix_Fields(mixing_type,field1,dfield1,iter_status)
!
use Iterate
use global
     implicit none

     integer, intent(in)          :: mixing_type
     !double precision, intent(in) :: lambda
     integer, intent(inout)       :: iter_status

     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field1
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: dfield1
   !  double precision, dimension(gridx,gridy,gridz,monomer_types):: dfield
!
     integer, parameter :: nvar = gridx*gridy*gridz*monomer_types
     double precision, save,dimension(nvar,0:mixing_dim)::var
     double precision, save, dimension(nvar,0:mixing_dim)::dvar

     integer :: x,y,z,n,monomer_type

if (iter_status==1)then
dvar(:,:) = 0.0d0
var=0.0d0
endif

   ! print*,sum( dvar(:,2)),iter_status,'up'
      n=0
      do monomer_type = 1,monomer_types
        do z = 1,gridz
         do y = 1,gridy
          do x = 1,gridx
            n=n+1
             var(n,0)  = field1(x,y,z,monomer_type)
             dvar(n,0) = dfield1(x,y,z,monomer_type)
          end do
         end do
        end do
      end do
      !  print*,( dvar(:,0)),iter_status,'m1'

      mixing: select case (mixing_type)
           case (1)
              call mixing_simple(var,dvar,nvar,lambda,iter_status)
           case (2)
              call mixing_lambda(var,dvar,nvar,lambda,iter_status)
           case (3)
              call mixing_anderson(var,dvar,nvar,lambda,iter_status)
                    
           case default
              print*, 'Invalid mixing option in find_saddle_point!'
              stop
      end select mixing
 !print*,sum( dvar(:,1)),iter_status,'down'
      n=0
      do monomer_type = 1,monomer_types
       do z = 1,gridz
         do y = 1,gridy
          do x = 1,gridx
            n=n+1
             field1(x,y,z,monomer_type) = var(n,0)
             dfield1(x,y,z,monomer_type) = dvar(n,0)
          end do
         end do
        end do
      end do
             !  print*,sum( dvar(:,:))!print*,sum( dvar(:,2)),iter_status
  ! print*,sum(dvar(:,1)),'up'
      iter_status = iter_status + 1

    return

  end subroutine Mix_Fields
  SUBROUTINE init_random_seed()

        implicit none
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)

    END SUBROUTINE   
    
    subroutine random_stdnormal(rand_F,sigma,mu)
   use Fourier_fftw3
   use global
   implicit none
   
   !standard deviation sigma
   !mean mu
   double precision,intent(in)::sigma,mu
   double complex, dimension(gridx/2+1,gridy, gridz,monomer_types,3), intent(out) :: rand_F
   double precision, dimension(gridx,gridy, gridz):: x_R
   integer::x,y,z,dimention,monomer_type
   !double precision,parameter :: pi=3.14159265
   double precision :: u1,u2

do monomer_type=1,monomer_types
   do dimention=1,3
      x_R=0
      do z = 1,gridz
         do y = 1,gridy
           do x = 1,gridx
   call random_number(u1)
   call random_number(u2)
  !u1=grnd()
  !u2=grnd()

   x_R(x,y,z) = dsqrt(-2*dlog(u1))*dcos(2*pi*u2)
   ! Transform to desired mean and standard deviation
   x_R(x,y,z) = mu + sigma * x_R(x,y,z)
   
   enddo
   enddo
   enddo

    call Real2Fourier(x_R(:,:,:),rand_F(:,:,:,monomer_type,dimention))
    enddo
    enddo
    return

end subroutine random_stdnormal


    subroutine random_stdnormal_real(x_r,sigma,mu,gridx,gridy,gridz)
   implicit none
   !standard deviation sigma
   !mean mu
   double precision,intent(in)::sigma,mu
   !double complex, dimension(gridx/2+1,gridy, gridz,3), intent(out) :: rand_F
   double precision, dimension(gridx,gridy, gridz),intent(out):: x_R
   integer::x,y,z,gridx,gridy,gridz
   !double precision,parameter :: pi=3.14159265
   double precision :: u1,u2
      x_R=0
  ! do k=1,3
do z = 1,gridz
       do y = 1,gridy
           do x = 1,gridx
   call random_number(u1)
   call random_number(u2)
   !u1=grnd() 
  ! u2=grnd() 
   x_R(x,y,z) = dsqrt(-2*dlog(u1))*dcos(2*pi*u2)
   ! Transform to desired mean and standard deviation
   x_R(x,y,z) = mu + sigma * x_R(x,y,z)
   
  ! enddo
enddo
   enddo
      enddo
    return
end subroutine random_stdnormal_real


!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     T H E    E N D 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end program main_blend
