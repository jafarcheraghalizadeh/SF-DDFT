!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM (MULTIBLOCK COPOLYMER MELT)
!
!       ALTERNATIVELY
!       
!          SCF CALCULATION OR
!          DYNAMIC DENSITY FUNCTIONAL CALCULATION WITH IMPOSED MOBILITY MATRIX
!
!          References (for DDFT): 
!              S. Mantha, S. Qi, F. Schmid, Macromolecules 53, 3409 (2020)
!                 doi:10.1021/acs.macromol.0c00130
!              F. Schmid, B. Li, Polymers 12, 2205 (2020)
!                 doi:10.3390/polym12102205
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_multiblock
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	    DECLARATIONS 
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    use Global
    use Molecule_Multiblock_Copolymer
    use Interactions_FloryHuggins
    use mtmod
    use Fourier_fftw3
     use Propagate
  
    implicit none
  
  ! GLOBAL PARAMETERS: DECLARED / DEFINED IN MODULE GLOBAL
  !
  ! double precision, parameter :: pi = 3.141592653589793d0
  ! integer, parameter  :: gridx, gridy, gridz 		! lattice points in each direction
  ! integer, parameter  :: ngrid                      ! total number of lattice points
  ! double precision    :: sizex, sizey, sizez     	! system size
  ! double precision    :: volume	               		! system volume
  ! double precision    :: dvol 				    ! volume element
  ! double precision    :: ds		                ! discretization along chain
  ! integer, paramter   :: mixing_dim             ! Mixing dimension for Anderson mixing
  !
    type(multiblock_copolymer) :: multiblock
  !
    double precision    ::  C      ! dimensionless polymer concentration
    double precision    ::  mu     ! chemical potential of polymer
  
    double precision    ::  energy
  
    integer             ::  itermax_scf    ! maximum number of iterations
    integer             ::  iterations     ! actual number of iterations
    double precision    ::  accuracy_scf   ! desired accuracy in iterations
    double precision    ::  accuracy       ! actual accuracy
   ! double precision    ::  lambda         ! mixing parameters
    integer             ::  mixing_type    ! type of mixing
  !      -------------  1: simple mixing, 2: lambda mixing, 3: anderson mixing ----------
  
    double precision    ::  wcutoff        ! cutoff for DDFT calculations (numerical stability)
  
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: density
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: field
  !
  ! REMOVE LATER
    double precision, dimension(gridx,gridy,gridz,monomer_types) ::  newfield, dfield
  ! REMOVE LATER
  
  !  double precision :: interaction_energy
  
   ! double precision, dimension(gridx,gridy,gridz,monomer_types) :: potential
  
    double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: mobility_qq
   ! double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: mobility
    double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: lambda_a
  
    double precision, dimension(gridx/2+1,gridy,gridz) :: q1
  
    logical           :: qexist, qread, qddft, qdmin
    character*24      :: file_name
    integer           :: polymer_label
    integer           :: polymer_ensemble   ! canonical:1, grand canonical:2
    integer           :: qinit   ! 1: cosine, 2: random densities, 3: random fields
  !  integer,parameter :: nblocks = 10 
  !   integer,parameter :: nblocks = 6
  !  integer,parameter :: nblocks = 7
    integer,parameter :: nblocks = 2
    double precision, dimension(nblocks) :: block_length
    integer, dimension(nblocks)          :: block_monomer_type
  !
    integer          :: number_of_timesteps ,start_time1,end_time,count_rate! number of time steps
    integer::start_time2,end_time2
    double precision :: timestep            ! length of time step
    double precision :: start_time          ! start time
  
    double precision :: dx, dy, dz, dum
    integer          :: io_error, x,y,z, step, i,i0, gcfstat,sample
    double precision::r
        double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier,dummyfourier1
        character(len=20) :: filenameq,filenameq_dens
  
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	     INITIALIZATIONS
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! ------ READ INPUT PARAMETERS 
   call system_clock(count_rate=count_rate)
  call system_clock(start_time1)
   call get_command_argument(1, filenameq)
   call get_command_argument(2, filenameq_dens)
   open(unit=10,file=trim(filenameq), status='replace',action='write',iostat=io_error)
   open(unit=100,file='E_'//trim(filenameq), status='replace',action='write',iostat=io_error)
  
   call init_random_seed()
   file_name='input_general'
   call read_input_parameters(file_name)
  !
  ! ------ GENERALvim INITIALIZATIONS 
  !
    call Initialize
  !
  ! ------ INITIALIZE MULTIBLOCK COPOLYMER
  !
      polymer_label = 1
      block_monomer_type = (/1,2/)
  !    block_monomer_type = (/1,2,1,2,1,2/)
  !    block_monomer_type = (/1,2,1,2,1,2,1/)
  !    block_monomer_type = (/1,2,1,2,1,2,1,2,1,2/)
       block_length = (/0.5d0,0.5d0/)
  !    block_length = (/0.5d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)
  !    block_length = (/0.1d0,0.1d0,0.1d0,0.4d0,0.1d0,0.1d0,0.1d0/)
  !    block_length = (/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)
      
      polymer_ensemble  = 1  ! canonical
     ! polymer_ensemble  = 2  ! grand canonical
      call Initialize_Multiblock(multiblock, polymer_label, nblocks, &
                       block_monomer_type, block_length, &
                       polymer_ensemble,C,mu)
  !
  ! ------ INITIALIZE INTERACTIONS
  !
    call Initialize_Interactions
  !
  ! ------ READ AND INITIALIZE MOBILITIES 
  !
    if (qddft) then
      !file_name='input_mobilities-KG400'
       file_name='Mobility-delta-print.txt'
      !file_name='input_mobilities'
    !  file_name= 'input_mobilities_d'
      call Read_Mobilities(file_name)
    end if
  !
  ! ------ CHECKS
  !
    call Calculate_Densities(field,density)
  
   call Check_Volume_Fraction(density)  ! for fully canonical with Flory Huggins interactions
  
  !
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
  !      INITIAL CONFIGURATION
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
    if (qread) then
       open (unit=33, file='cnf.in', status='old', action='read', iostat=io_error)
       if(io_error /=0) then
         write(*,*) 'Error', io_error , ' while trying to open cnf.in'
       end if
       read (33,*) field
       close (unit=33)
       call Calculate_Densities(field,density)
     else
      if (qinit.eq.1) then
        do x = 1,gridx
           do z = 1,gridz
             density(x,:,z,1) = multiblock%rho &
  !            * ( block_length(1)+block_length(3)+block_length(5) ) &
               *  block_length(1) &
  !             * ( 1.0d0 + 0.2d0*cos(dble(x)/dble(gridx)*2*pi) )
               * ( 1.0d0 - 0.2d0*cos(dble(x)/dble(gridx)*2*pi)*cos(dble(z)/dble(gridz)*2*pi) )
             density(x,:,z,2) = 1.0d0 - density(x,:,z,1)
           end do
         end do
  !       call Get_Conjugate_Fields(density,field,1.d-6,100000,3,gcfstat)
         call Calculate_Mean_Fields(density,field)
         call Calculate_Densities(field,density)
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
  
          density(:,:,:,1) = multiblock%rho &
  !           *(block_length(1)+block_length(3)+block_length(5)) &
             *  block_length(1) &
             *(1.0d0 + (2.0d0*density(:,:,:,1) -1.0d0)*0.01d0)
          density(:,:,:,2) = 1.0d0 - density(:,:,:,1)
  !             --- densities not yet normalized
         call Get_Conjugate_Fields(density,field,0.01d0,10000,3,gcfstat)
         call Calculate_Densities(field,density)
         call Get_Conjugate_Fields(density,field,0.01d0,10000,3,gcfstat)
   
  !             --- continue with normalized densities
       else if (qinit .eq. 3) then
          do i = 1,monomer_types
            do z=1,gridz
              do y=1,gridy
                 do x=1,gridx
                 call random_number(r)
                   field(x,y,z,i) =r!grnd()
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
                  end do
              end do
            end do
          print*, 'Successful reading of density'
          close (unit=165)
    call Get_Conjugate_Fields(density,field,0.01d0,10000,3,gcfstat)
          call Calculate_Densities(field,density)
          call Get_Conjugate_Fields(density,field,0.01d0,10000,3,gcfstat)
       else
         print*, 'Initialization not well-defined!'
         stop
       end if
     end if
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	   OUTPUT OF INITIAL CONFIGURATION
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
     !   open(unit=60,file='initial_dens.dat', status='replace',action='write',iostat=io_error)
       ! open(unit=61,file='initial_field.dat', status='replace',action='write',iostat=io_error)
          !dx = sizex/float(gridx)
         ! dy = sizey/float(gridy)
         ! dz = sizez/float(gridz)
         ! do x=1,gridx
             !do y=1,gridy
               ! do z=1,gridz
                 !  write(60, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                        !       (density(x,y,z,i),i=1,monomer_types)
               !    write(61, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                      !         (field(x,y,z,i),i=1,monomer_types)
               ! end do
            ! end do
            ! write(60, *)
            ! write(61, *)
         ! end do
      !  close(unit=60)
       ! close(unit=61)
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !      DDFT SIMULATION / SCF CALCULATION
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    if (qddft) then
  
      print*, 'BEGIN DDFT SIMULATION'
      write(30,*) 'BEGIN DDFT SIMULATION'
      write(30,*)
  
      qdmin = .true.
  
      do step = 1,number_of_timesteps
        call DDFT_Timestep_Euler(density,field,lambda_a,timestep,qdmin)
  
        if (((step/50)*50).eq. step) then                 
         call Calculate_Energy    
         write(30,*) step, energy ! FOR TESTING
         write(*,'(F10.3,I10,A,F10.4,F10.4,F10.4)')energy, step/50, ', t=',timestep*step
        end if                                                
  !      
        if (((step/2000)*2000).eq. step) then                 ! FOR TESTING
         call Calculate_Energy                                ! FOR TESTING
       !!!!  write(30,*) step, energy, multiblock%rho, multiblock%mu  ! FOR TESTING
        ! call Check_Volume_Fraction(density)                  ! FOR TESTING
        end if                                                ! FOR TESTING
  !
        i = nint(start_time/timestep) + step
        i0 = nint(0.1/timestep)
        if (((step/50)*50).eq. step .and. step>1) then                   ! FOR TESTING
          dx = sizex/float(gridx)
          dy = sizey/float(gridy)
          dz = sizez/float(gridz)
          do z=1,gridz
             do y=1,gridy
                do x=1,gridx
                write(10,*)density(x,y,z,1),density(x,y,z,2),field(x,y,z,1),field(x,y,z,2)
                end do
             end do
          end do
                 write(100,*)timestep*step,energy
                 
         ! open(unit=33, file='cnf.out', status = 'replace', action='write', iostat=io_error)
         !    write (33,*) field
         ! close (unit=33)
  
        end if
  
  
      end do
   else
  
      print*, 'BEGIN SCF CALCULATION'
      write(30,*) 'BEGIN SCF CALCULATION'
      write(30,*)
  
      call Find_SCF_Saddle_Point(density,field)
      call Calculate_Densities(field,density)
  
    end if
  
    call Calculate_Energy
  
    if (qddft) then
      write(30,*) '   Final values after ', number_of_timesteps, 'steps'
      print*, ' Final values after ', number_of_timesteps, 'steps'
    else
      write(30,*) '   Final values '
      print*, ' Final values '
  
    end if
    write(30,*) '   Energy: ', energy, ' C: ', multiblock%rho, ' mu_scf: ', multiblock%mu
    write(30,*)
  !
    print*, '   Energy: ', energy, energy/volume, ' C: ', multiblock%rho, ' mu_scf: ', multiblock%mu
  !
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	   OUTPUT OF FINAL CONFIGURATION
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
       ! open(unit=60,file='final_dens.dat', status='replace',action='write',iostat=io_error)
       ! open(unit=61,file='final_field.dat', status='replace',action='write',iostat=io_error)
        !  write(60,*) "# Energy : ", energy
          !write(60,*) 
          !write(61,*) "# Energy : ", energy
          !write(61,*) 
         ! dx = sizex/float(gridx)
         ! dy = sizey/float(gridy)
         ! dz = sizez/float(gridz)
        !  do x=1,gridx
           !  do y=1,gridy
               ! do z=1,gridz
                  ! write(60, *)   x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                            !   (density(x,y,z,i),i=1,monomer_types)
                   !write(61, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                             !  (field(x,y,z,i),i=1,monomer_types)
               ! end do
          !   end do
             !write(60, *)
           !  write(61, *)
          !end do
       ! close(unit=60)
       ! close(unit=61)
  !
     !   open(unit=33, file='cnf.out', status = 'replace', action='write', iostat=io_error)
         !  write (33,*) field
     !   close (unit=33)
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	   FINISH (CLOSE FILES ETC.)
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   call system_clock(end_time)
   print*, 'Elapsed time: ', real(end_time - start_time1)/count_rate
    close (unit=30)
    call finish
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	   END OF MAIN PROGRAM. NOW INNER SUBROUTINES
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !end do
  contains
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	   SUBROUTINE TO READ INPUT PARAMETERS  ! THIS DEPENDS ON SYSTEM !!!
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine Read_Input_Parameters(file_name)
      
      implicit none
  
      character*24, intent(in) :: file_name
      integer                  :: io_error
  
      open(unit=20,file=file_name,status='old',action='read',iostat=io_error)
  
        if(io_error ==0) then
  
          read(20,*) sizex
          read(20,*) sizey
          read(20,*) sizez
          read(20,*)
          read(20,*) C
          read(20,*) mu
          read(20,*)
          read(20,*) ds
          read(20,*)
          read(20,*) start_time  ! start time
          read(20,*) number_of_timesteps
          read(20,*) timestep
          read(20,*)
          read(20,*) qread  ! qread = .true. : read from file cnf.in
          read(20,*) qinit  ! only for qread=.false.  
              ! 1=cosine, 2=random densities, 3=random fields
          read(20,*) qddft  ! qepd = .true. : EPD simulation; otherweise SCF calculation
          read(20,*)
          read(20,*) wcutoff
  
        else
          write(*,*) 'Error', io_error , ' while trying to open ', file_name
        end if
  !
      close(unit=20)
  
    end subroutine Read_Input_Parameters
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	   SUBROUTINE TO READ AND INITIALIZE MOBILITIES  ! THIS DEPENDS ON SYSTEM !!!
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine Read_Mobilities(file_name)
      
      use Fourier_fftw3
  
      implicit none
  
      character*24, intent(in) :: file_name
  
      integer, parameter       :: nqmax=1000
      double precision, dimension(nqmax,monomer_types,monomer_types) :: mobility_table
      double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: mobility
      
      double precision,allocatable::q_values(:)
      double precision :: dum, dq, qq, q, small
      integer          :: io_error, nq, iq, x, y, z, iqcut,type1,type2,idx
      real:: a,b,c
      double precision::q_low,q_high,d_low,d_high
  !   ------------ read mobilitiesMobility-delta-print.txt 
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
             read(20,*) dum , mobility_table(iq,1,1) &
                       ,mobility_table(iq,1,2),mobility_table(iq,2,2)
           !  read(20,*) dum , mobility_table(iq,1,1),mobility_table(iq,1,2)
              q_values(iq) = dum
            ! if (abs(dum -iq*dq) .ge. small ) then
              ! print*, 'q values in mobility are not consistent', dum, iq*dq
              ! stop 
            ! end if
          end do
          mobility_table(:,2,1) = mobility_table(:,1,2)
        !  mobility_table(:,2,2) = mobility_table(:,1,1)
        ! mobility_table=mobility_table*0.8d0
        else
          write(*,*) 'Error', io_error , ' while trying to open ', file_name
        end if
      close(unit=20)
  mobility=0.0d0
  mobility_qq=0.0d0
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
                      do while (idx <= nq .and. q_values(idx) <= q)
  
                                idx = idx + 1
                                if (idx>nq) exit
                      end do
                      idx = idx - 1
                  if (idx <= 1) then
          mobility(x,y,z,type2,type1) = mobility_table(1,type2,type1)
                    elseif (idx >= nq) then
                        mobility(x,y,z,type2,type1) = mobility_table(nq,type2,type1)
                    else
                        q_low  = q_values(idx)
                        q_high = q_values(idx + 1)
                        d_low  = mobility_table(idx,type2,type1)
                       d_high = mobility_table(idx + 1,type2,type1)
                      
                  ! Perform linear interpolation
                        mobility(x,y,z,type2,type1) = d_low + (q - q_low) * (d_high - d_low) / (q_high - q_low)
                       !  print*, mobility(x,y,z,type2,type1)
        if ((q_high - q_low)==0)print*,'ERROR Mobility: Zero in denominator'
                    end if
              
              
   if ((mobility(x,y,z,type2,type1)>=0))then
         mobility_qq(x,y,z,type2,type1) = min(qq*mobility(x,y,z,type2,type1),wcutoff)
          else if ((mobility(x,y,z,type1,type2)<0))then
         mobility_qq(x,y,z,type2,type1) = max(qq*mobility(x,y,z,type2,type1),-wcutoff+1.0e-30)
         endif
  
            
  
        end do
       end do
      end do
      enddo
      enddo
      lambda_a(:,:,:,:,:) = 0.0d0
  
      ! Calculate lambda_a(:,:,:,1,1)
   
  
      ! Apply conditions where lambda_a(:,:,:,1,1) is not zero
      do x = 1, gridx/2+1
          do y = 1, gridy
              do z = 1, gridz
                   lambda_a(x,y,z,1,1) = dsqrt(mobility_qq(x,y,z,1,1))
                  if (lambda_a(x,y,z,1,1) > 0.0d0 ) then
                      ! Calculate lambda_a(:,:,:,2,1)
                    !
                      !print*,mobility_qq(x,y,z,2,1), (lambda_a(x,y,z,1,1)),x,y,z
                      lambda_a(x,y,z,2,1) = mobility_qq(x,y,z,2,1) / (lambda_a(x,y,z,1,1))
                  else
                  continue
                  end if
                      ! Calculate lambda_a(:,:,:,2,2)
                      lambda_a(x,y,z,2,2) = dsqrt(mobility_qq(x,y,z,2,2) - lambda_a(x,y,z,2,1)**2)
  
                      write(150,*) mobility(x,y,z,1,1)
              end do
          end do
      end do
  call flush(150)
     !=====================================================================================================
    end subroutine Read_Mobilities
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
  !
      implicit none
  
      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout)  :: field
      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: density
  
      double precision, dimension(gridx,gridy,gridz,monomer_types) :: expfield
  
      integer :: monomer_type,i
  
      call Calculate_Expfield(field,expfield)
      call Propagate_Multiblock(multiblock,expfield)
  
      density = 0.0d0
      do i = 1,nblocks
        monomer_type = multiblock%monomer_type(i)
        density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + multiblock%density(:,:,:,i)   
      end do
      
  
  !    density(:,:,:,diblock%monomer_type) = diblock%density
  
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
      double precision              :: chain_energy
  
      chain_energy = multiblock%energy - dvol*sum(density*field)
      energy = chain_energy + interaction_energy(density)
  
      return
  
    end subroutine Calculate_Energy
  !
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !      SUBROUTINE TO CALCULATE FIELDS CORRESPONDING TO GIVEN DENSITIES  
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine Get_Conjugate_Fields(density,field,accuracy_scf,itermax_scf,mixing_type,gcfstat)
  !
     use Propagate
        use, intrinsic :: ieee_arithmetic
       implicit none
  
       double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
       double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)    :: density
       integer, intent(out) :: gcfstat
   
       double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
       double precision, dimension(gridx,gridy,gridz,monomer_types) ::  newdensity
  
       double precision, intent(in) :: accuracy_scf
       integer, intent(in) :: itermax_scf, mixing_type
   
       double precision :: accuracy, small
       integer :: iter_status, iter
   
  !    accuracy_scf = 0.001
  !    itermax_scf = 100000
  !    mixing_type = 3
  
       
  
       iter_status = 1
   
       do iter = 1,itermax_scf
  
   ! if (norm(2)-norm(1)>0 .and. iter>1)then
     !   lambda=max(lambda_min,0.9*lambda)
    !else
      !  lambda=min(lambda_max,1.1*lambda)
   ! end if
  lambda=0.40d0
  
         call Calculate_Densities(field,newdensity)
         dfield = newdensity - density
         call Check_Accuracy(dfield,accuracy)
  
         if (.not.(accuracy > accuracy_scf*0.5) ) then  ! double negation to deal with NaNs
         ! PRINT*,'exit Get_Conjugate_Field after', iter, 'steps '
         ! PRINT*,'Accuracy: ', accuracy
           exit
         end if
         iter_status = min(iter,mixing_dim+1)
         
         
         call Mix_Fields(mixing_type,field,dfield,iter_status)
  
       end do
  
       if (.not.(accuracy < accuracy_scf*0.5)) then  ! double negation to deal with NaNs
         if (.not.(accuracy <  accuracy_scf)) then  ! double negation to deal with NaNs
          ! print*, 'Warning: Conjugate fields not found within ',iter,' iteration steps!'
         !  print*, 'Accuracy: ', accuracy
           gcfstat = 2 ! iteration not successful
           if (.not.(accuracy <  200.0d0)) then  ! double negation to deal with NaNs
             print*, "accuracy > 200.0"
             stop
           end if
         else
          ! print*, 'Caution: Finding conjugate fields takes long!'
         !  print*, 'Accuracy after ', iter, 'iteration steps: ', accuracy
          ! print*,maxval(density),maxval(newdensity)
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
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !      SUBROUTINE TO FIND SCF SADDLE POINT
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
     subroutine Find_SCF_Saddle_Point(density,field)
  !
      use Propagate
  !   use Interactions_FloryHuggins
  
       implicit none
  
       double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
       double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density
   
       double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
       double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity, newfield
   
       double precision ::  accuracy, accuracy_scf
       integer :: mixing_type, iter_status, iter, iter_max
   
       accuracy_scf = 0.00001
       mixing_type = 3
       iter_max = 1000000
     !  lambda = 0.2
   
       iter_status = 1
   
       do iter = 1,iter_max
   
         call Calculate_Densities(field,newdensity)
         call Calculate_Mean_Fields(newdensity,newfield)
         dfield = newfield - field
         call Check_Accuracy(dfield,accuracy)
        if ( ( (iter/400)*400) .eq. iter) then
          print*, iter, accuracy
        end if
   
         if (.not.(accuracy > accuracy_scf) ) then  ! double negation to deal with NaNs
          ! print*, 'SCF saddle point found after ',iter,' iteration steps!'
         !  print*, 'Accuracy: ', accuracy
           exit
         end if
         iter_status = min(iter,1000)
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
  !	   SUBROUTINE TO MIX FIELDS IN SCF ITERATION
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine Mix_Fields(mixing_type,field,dfield,iter_status)
  !
     use Iterate
  !
       implicit none
  
       integer, intent(in)               :: mixing_type, iter_status
      ! double precision, intent(in)      :: lambda
        
       double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
       double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: dfield
  !
       integer, parameter :: nvar = gridx*gridy*gridz*monomer_types
       double precision, dimension(nvar,0:mixing_dim), save :: var,dvar 
  
       integer :: x,y,z,n,monomer_type
  
        n=0
        do monomer_type = 1,monomer_types 
          do z = 1,gridz
           do y = 1,gridy
            do x = 1,gridx
              n=n+1
               var(n,0)  = field(x,y,z,monomer_type)
               dvar(n,0) = dfield(x,y,z,monomer_type)
            end do
           end do
          end do
        end do
  
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
  
        n=0
        do monomer_type = 1,monomer_types 
          do z = 1,gridz
           do y = 1,gridy
            do x = 1,gridx
              n=n+1
               field(x,y,z,monomer_type) = var(n,0)  
               dfield(x,y,z,monomer_type) = dvar(n,0) 
            end do
           end do
          end do
        end do
  
      return
  
    end subroutine Mix_Fields
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !           CALCULATE ACHIEVED ACCURACY
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine Check_Accuracy(dfield,accuracy)
  
      integer :: x,y,z
      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: dfield
      double precision, intent(out)       :: accuracy
  
     ! accuracy=sum(Abs(dfield))*dvol
      accuracy=maxval(Abs(dfield))
  !if( maxval(Abs(dfield))>2.0d0)then
  !print*,'accuracy:',accuracy
  
  !end if 
      return
    end subroutine check_accuracy
  !
  !
  !
  
  
  
  
  
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
  
        double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density
        double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
  
        !double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types), intent(in) :: mobility_qq
        !double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types), intent(in) :: mobility
        double precision, intent(in) :: timestep
        logical, intent(in)          :: qdmin
  
        double precision, dimension(gridx,gridy,gridz,monomer_types)   :: newfield
        double precision, dimension(gridx,gridy,gridz,monomer_types)   :: mulocal
        double precision, dimension(gridx,gridy,gridz,monomer_types)   :: ddensity
  
        double complex, dimension(gridx/2+1,gridy,gridz,monomer_types) :: mu_q
        double complex, dimension(gridx/2+1,gridy,gridz,monomer_types) :: ddensity_q
        double complex, dimension(gridx/2+1,gridy,gridz,monomer_types) :: j_a
        double precision, dimension(gridx,gridy,gridz,monomer_types)   :: j_a_r
        double precision, dimension(gridx,gridy,gridz)   :: dummy
        double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier
        double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types),intent(in) :: lambda_a
  !
        double precision :: dsmall
        integer :: x,y,z,n,monomer_type,type1,type2, gcfstat
        double complex, dimension(gridx/2+1,gridy, gridz,monomer_types) :: rand_F
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
  !
        do monomer_type = 1, monomer_types
          dummy = mulocal(:,:,:,monomer_type)
          call Real2Fourier(dummy,dummyfourier)
  !        dummyfourier = - laplace*dummyfourier
          mu_q(:,:,:,monomer_type) = dummyfourier
        end do
     !  do z = 1,gridz
             !   do y = 1,gridy
           ! do x = 1,gridx/2+1
           ! write(250,*)  real(mu_q(x,y,z,1)), aimag(mu_q(x,y,z,1)),&
           ! real(mu_q(x,y,z,2)), aimag(mu_q(x,y,z,2))
           ! enddo
           ! enddo
           ! enddo
  !     Density derivative in Fourier space
  !
       
       ddensity_q = 0.0d0
       j_a= 0.0d0
     !  lambda_a(:,:,:,1,1)=sqrt(mobility(:,:,:,1,1))
      ! lambda_a(:,:,:,2,1)=mobility(:,:,:,2,1)/lambda_a(:,:,:,1,1)
      ! lambda_a(:,:,:,2,2)=sqrt(mobility(:,:,:,2,2)- lambda_a(:,:,:,2,1)**2)
  
       do type1=1,monomer_types
            do type2=1,monomer_types
                    ddensity_q(:,:,:,type1) =  ddensity_q(:,:,:,type1) &
                     -mobility_qq(:,:,:,type2,type1)*mu_q(:,:,:,type2)
             enddo
              enddo
         !    ddensity_q(:,:,:,1) =  ddensity_q(:,:,:,1)-mu_q(:,:,:,1)
           !  ddensity_q(:,:,:,2) =  ddensity_q(:,:,:,2)-mu_q(:,:,:,2)
  
     call random_stdnormal_FFT(rand_F,kt,0.0d0)
     !call random_stdnormal_FFT(rand_F1,kt,0.0d0,gridx,gridy,gridz)
         ! call random_stdnormal_real(ran, 1.0d0,0.0d0,gridx,gridy,gridz)
  
    ! J_1=============================================
       !	j_a(x,y,z,1)=j_a(x,y,z,1)&
         !+(rand_F(x,y,z,1))*nabla(x,y,z,1)* lambda_a(x,y,z,1,1)+(rand_F(x,y,z,2))&
         !*nabla(x,y,z,2)* lambda_a(x,y,z,1,1)+(rand_F(x,y,z,3))*nabla(x,y,z,3)*lambda_a(x,y,z,1,1)
      j_a(:,:,:,1)=j_a(:,:,:,1)+rand_F(:,:,:,1)* lambda_a(:,:,:,1,1)
         
          ! J_2============================================
   !	j_a(x,y,z,2)=j_a(x,y,z,2)&
     !+((rand_F(x,y,z,1))*nabla(x,y,z,1)* lambda_a(x,y,z,2,1)+(rand_F(x,y,z,2))&
     !*nabla(x,y,z,2)* lambda_a(x,y,z,2,1)+(rand_F(x,y,z,3))*nabla(x,y,z,3)* lambda_a(x,y,z,2,1))
     j_a(:,:,:,2)=j_a(:,:,:,2)+rand_F(:,:,:,1)* lambda_a(:,:,:,2,1)
    
         !j_a(x,y,z,2)=j_a(x,y,z,2)&
         !+((rand_F1(x,y,z,1))*nabla(x,y,z,1)* lambda_a(x,y,z,2,2)+(rand_F1(x,y,z,2))&
         !*nabla(x,y,z,2)* lambda_a(x,y,z,2,2)+(rand_F1(x,y,z,3))*nabla(x,y,z,3)*lambda_a(x,y,z,2,2))
           
      j_a(:,:,:,2)=j_a(:,:,:,2)+rand_F(:,:,:,2)* lambda_a(:,:,:,2,2)
  
  !     Density derivative in Real space
  
        do monomer_type = 1, monomer_types
          dummyfourier = ddensity_q(:,:,:,monomer_type)
          call Fourier2Real(dummy,dummyfourier)
          ddensity(:,:,:,monomer_type) = dummy
         
          dummyfourier = j_a(:,:,:,monomer_type)
  
          call Fourier2Real(dummy,dummyfourier)
          j_a_r(:,:,:,monomer_type)=dummy*sqrt(2.0d0)*sqrt(timestep)
  
  !        offset(monomer_type) = sum(ddensity(:,:,:,monomer_type))/ngrid
  !        ddensity(:,:,:,monomer_type) = ddensity(:,:,:,monomer_type) &
  !            + offset(monomer_type)
        end do
  
  !     Euler step
        dsmall = 1.d-4
       do monomer_type = 1, monomer_types
        density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + &
        timestep*(ddensity(:,:,:,monomer_type)) + j_a_r(:,:,:,monomer_type)
        enddo
        
        call Get_Conjugate_Fields(density,field,0.01d0,2000,3,gcfstat)
     
        if (qdmin .and. (gcfstat.ge.1)) then
         call Calculate_Densities(field,density) 
        end if
  !
       return
   
     end subroutine DDFT_Timestep_Euler
  
  
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
      
      subroutine random_stdnormal_FFT(rand_F,sigma,mu)
      use Fourier_fftw3
      use global
      implicit none
     !standard deviation sigma
     !mean mu
     double precision,intent(in)::sigma,mu
     double complex, dimension(gridx/2+1,gridy, gridz,monomer_types), intent(out) :: rand_F
     double precision, dimension(gridx,gridy, gridz,monomer_types):: x_R
     integer::x,y,z,monomer_type
     !double precision,parameter :: pi=3.14159265
     double precision :: u1,u2
  
     x_R=0
     do monomer_type = 1, monomer_types
      do z = 1,gridz
         do y = 1,gridy
        do x = 1,gridx
     call random_number(u1)
     call random_number(u2)
     !u1=grnd()
     !u2=grnd()
     x_R(x,y,z,monomer_type) = dsqrt(-2*dlog(u1))*dcos(2*pi*u2)
     ! Transform to desired mean and standard deviation
     x_R(x,y,z,monomer_type) = mu + sigma * x_R(x,y,z,monomer_type)
        enddo
        enddo
     enddo
     
      call Real2Fourier((x_R(:,:,:,monomer_type)),rand_F(:,:,:,monomer_type))
  
  enddo
      return
  end subroutine random_stdnormal_FFT
  
  !
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !	     T H E    E N D 
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  
  end program main_multiblock
  