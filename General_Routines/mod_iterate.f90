!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR ITERATION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Iterate

  use Global
  implicit none

  !integer, parameter :: mixing_dim = 2

contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           DIFFERENT MIXING ROUTINES
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!		1. SIMPLE MIXING
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mixing_Simple(field,dfield,fieldsize,lambda,iter)

    double precision, intent(in) :: lambda
    integer, intent(in)          :: fieldsize
    integer, intent(in)          :: iter  ! not used in this routine
    double precision, dimension(fieldsize,0:mixing_dim), intent(inout) :: field
    double precision, dimension(fieldsize,0:mixing_dim), intent(inout) :: dfield
  
    field(:,0) = field(:,0) + lambda*dfield(:,0)

    return
  end subroutine Mixing_Simple
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!		2. MIXING BY VARIABLE LAMBDA (FOR ITER > 1)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mixing_Lambda(field,dfield,fieldsize,lammax,iter)

    double precision, intent(in) :: lammax
    integer, intent(in)          :: fieldsize
    integer, intent(in)          :: iter
    double precision, dimension(fieldsize,0:mixing_dim), intent(inout) :: field
    double precision, dimension(fieldsize,0:mixing_dim), intent(inout) :: dfield

    double precision, dimension(fieldsize)                :: ddiff

    double precision, save  :: lammin=0.00001d0

    double precision        :: lambda,lambda_num,lambda_nom
  
    if (iter == 1) then
      lambda = lammin
    else
      ddiff  = dfield(:,0) - dfield(:,1)
      lambda_num = sum( ( field(:,0)-field(:,1))**2 )
      lambda_nom = sum( (dfield(:,0)-dfield(:,1))**2 )
      if (lambda_nom > 1.0d-50) then
        lambda=sqrt(lambda_num/lambda_nom)
        if(lambda<lammin) then
          lambda=lammin
        else if (lambda>lammax) then
          lambda=lammax
      end if
      else
        lambda=lammin
      end if
    end if

    field(:,1)  = field(:,0)
    dfield(:,1) = dfield(:,0)

    field(:,0) = field(:,0) + lambda*dfield(:,0)

    return
  end subroutine Mixing_Lambda
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!		2. ANDERSON MIXING (FOR ITER > DIM)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mixing_Anderson(field,dfield,fieldsize,lam_mix,iter)
      use, intrinsic :: ieee_arithmetic
    double precision, save       :: lam00 = 0.001d0

    double precision, intent(in) :: lam_mix
    integer, intent(in)          :: fieldsize
    integer, intent(in)          :: iter
    double precision, dimension(fieldsize,0:mixing_dim), intent(inout) :: field
    double precision, dimension(fieldsize,0:mixing_dim), intent(inout) :: dfield

    double precision, dimension(fieldsize,mixing_dim)   :: ddiff
    double precision, dimension(mixing_dim,mixing_dim)  :: u_matrix
    double precision, dimension(mixing_dim)             :: v_vector, x_vector
    double precision, dimension(fieldsize)    ::   field_mix
    double precision, dimension(fieldsize)    ::   dfield_mix
    double precision :: lambda,hu
    integer :: i,j
 integer::x,y,z
    if (iter > mixing_dim) then

      do i = 1,mixing_dim
         ddiff(:,i) = dfield(:,0) - dfield(:,i)
  
      end do

      do i = 1,mixing_dim
         v_vector(i) = sum( ddiff(:,i)*dfield(:,0) )
         do j = 1,mixing_dim
            u_matrix(i,j) = sum( ddiff(:,i)*ddiff(:,j) )
         !  if (abs(u_matrix(i,j)) > huge(hu)) print*, 'umatrix is infinity'
         if (sum( dfield(:,0) )<-5 )print*, sum( dfield(:,0) ),sum( dfield(:,1) ),sum( dfield(:,2) )
        if (u_matrix(i,j)==0) print*, sum( dfield(:,0) ),sum( dfield(:,1) ),sum( dfield(:,2) )
         end do
      end do
           do x=1,gridx*gridy*gridz

       if (.not. ieee_is_finite(field(x,0)))then
       print*,'inf after'
       !stop
       end if
       end do
      call Solve_LinearSystem_LU(u_matrix,v_vector,x_vector,mixing_dim)

      field_mix  =  field(:,0) 
      dfield_mix = dfield(:,0)
      do i = 1,mixing_dim
         field_mix  =  field_mix - x_vector(i) * (field(:,0)-field(:,i))
         dfield_mix = dfield_mix - x_vector(i) * ddiff(:,i)
      end do

      lambda = lam_mix

    else

      field_mix = field(:,0)
      dfield_mix = dfield(:,0)

      lambda = lam00

    end if

    do i = mixing_dim,1,-1
       field(:,i)  =  field(:,i-1)
       dfield(:,i) = dfield(:,i-1)

    end do
    norm(2)=norm(1)
    field(:,0) = field_mix + lambda * dfield_mix
    
    !norm(1) = sum(abs( dfield_mix))
     norm(1) = sum(abs(field(:,0)- field_mix)) ! Measure the L1 norm of the field's new state
    
    
  end subroutine Mixing_Anderson
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINE TO SOLVE LINEAR EQUATION SYSTEM AX=B VIA CHOLESKY DECOMPOSITION
!	    A MUST BE SYMMETRIC AND POSITIVE DEFINITE
!           L IS A LOWER TRIANGULAR MATRIX WITH A = L^T L
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Solve_LinearSystem_Cholesky(A,B,X,d)

    integer, intent(in) :: d
    double precision, dimension(d,d), intent(in) :: A
    double precision, dimension(d),   intent(in) :: B
    double precision, dimension(d),  intent(out) :: X
    double precision, dimension(d,d) :: L
    double precision, dimension(d)     :: Y
    integer :: i,j,k
!
! ---------- CHOLESKY DECOMPOSITION OF MATRIX A
!
    L = 0.d0
!
    L(1,1) = sqrt(A(1,1))
    do i = 2,d
       do j = 1,i-1
          L(i,j) = ( A(i,j) - sum( L(i,1:(j-1))*L(j,1:(j-1)) ) ) / L(j,j)
       end do
       L(i,i) = sqrt( A(i,i) - sum( L(i,1:(i-1))**2 ) )
    end do
!
! ---------- SOLVE LY=B
!
    do i = 1,d
       Y(i) = ( B(i) - sum( L(i,1:(i-1))*Y(1:(i-1)) ) ) / L(i,i)
    end do
!
! ---------- SOLVE L^T X = Y
!
    do j = d,1,-1
       X(j) = ( Y(j) - sum( L((j+1):d,j)*X((j+1):d) ) ) / L(j,j)
    end do

    return
  end subroutine Solve_LinearSystem_Cholesky
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINE TO SOLVE LINEAR EQUATION SYSTEM AX=B VIA LU DECOMPOSITION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Solve_LinearSystem_LU(A,B,X,n)

    integer, intent(in) :: n
    double precision, dimension(n,n), intent(in) :: A
    double precision, dimension(n),   intent(in) :: B
    double precision, dimension(n),  intent(out) :: X

    integer, dimension(n)            :: indx
    double precision, dimension(n,n) :: aa
    double precision                 :: xx
!
! ---------- LU DECOMPOSITION OF MATRIX A
!
    X = B
    aa = A
!
    call ludcmp(aa,n,indx,xx)
    call lubksb(aa,n,indx,X)
!
    return
!
  end subroutine Solve_LinearSystem_LU
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           NUMERICAL RECIPE ROUTINES FOR LU DECOMPOSITION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     SUBROUTINE ludcmp(a,n,indx,d)
!
      integer, intent(in)                     :: n
      integer, dimension(n), intent(out)      :: indx
      double precision, intent(out)           :: d
      double precision, dimension(n,n),intent(inout) :: a
!
      integer, parameter          :: NMAX=500
      double precision, parameter :: tiny_number = 1.0d-20
!
      integer          :: i,imax,j,k
      double precision :: aamax,dum,dsum,vv(NMAX)
!
      d=1.
      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        end do
        if (aamax.eq.0.) then
        print*, 'singular matrix in ludcmp'
        print*,a
        end if
        vv(i)=1./aamax
      end do
!
      do j=1,n
!
        do i=1,j-1
          dsum=a(i,j)
          do k=1,i-1
            dsum=dsum-a(i,k)*a(k,j)
          end do
          a(i,j)=dsum
        end do
        aamax=0.
        do i=j,n
          dsum=a(i,j)
          do k=1,j-1
            dsum=dsum-a(i,k)*a(k,j)
          end do
          a(i,j)=dsum
          dum=vv(i)*abs(dsum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        end do
        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          end do
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=tiny_number
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          end do
        endif
      end do
!
      return
!
      END subroutine ludcmp
!
!-------------------------------------------------------
!
      SUBROUTINE lubksb(a,n,indx,b)
!
      integer, intent(in)                    :: n
      integer, dimension(n), intent(in)      :: indx
      double precision, dimension(n), intent(inout) :: b
      double precision, dimension(n,n), intent(in)  :: a

      integer          :: i,ii,j,ll
      double precision :: dsum
!
      ii=0
      do  i=1,n
        ll=indx(i)
        dsum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            dsum=dsum-a(i,j)*b(j)
          end do
        else if (dsum.ne.0.) then
          ii=i
        endif
        b(i)=dsum
      end do
!
      do i=n,1,-1
        dsum=b(i)
        do j=i+1,n
          dsum=dsum-a(i,j)*b(j)
        end do
        b(i)=dsum/a(i,i)
      end do
!
      return
!
      END subroutine lubksb
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module iterate
