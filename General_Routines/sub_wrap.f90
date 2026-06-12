!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    INITIALIZATIONS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize

    use Global
    use mtmod
    use Fourier_fftw3            
    use Propagate
!
    implicit none
!
! ------ FFTW INITIALIZATION -----------
!
    call Initialize_Fourier
!
! -------- TEST OF FOURIER TRANSFORM
!
!   call fourier_test
!
! ------- INITIALIZATION OF RNG --------
!
    call sgrnd(2007)
!  
! -------- LAPLACE OPERATOR ------------
!
    call Calculate_Laplace           
!
! ------ OTHER INITIALIZATIONS ---------
!
    call Calculate_Volume_Parameters
!
!--------------------------------------
    return

  end subroutine initialize
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    CLOSING ROUTINES
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Finish

    use Global
    use Fourier_fftw3            
!
    implicit none
!
! ------ FFTW CLOSING -----------
!
    call End_Fourier

  end subroutine finish
