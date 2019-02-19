!dir$ if .not. defined(noMKL)
!recursive communication interface for iterative sparse solver
include 'mkl_rci.f90'
!dir$ end if
    
module laWrapperLib
use constants
use arrayOpsLib
!dir$ if .not. defined(noMKL)
use lapack95
use mkl_rci
!dir$ end if
implicit none

    private
    
    !have emergency solver
    public:: solveGeneralLES
    public:: solveTridiagonalLES
    public:: triFactorSquareMat
    !--
    public:: polyfit
    public:: csrAx
    public:: fgmres
    public:: fgmresILUT
    
    !rely on the lapack
    !dir$ if .not. defined(noMKL)
    public:: choleskyFactorSym
    public:: solveLinearLeastSquare
    public:: solveSymmetryLES
    public:: eigenSymTriDiagonal
    public:: inverseGeneralSquareMat
    public:: eigenTriDiagonal
    !dir$ end if
    
    
!--------------------------------------------------------
    interface solveLinearLeastSquare
        procedure:: solveLinearLeastSquare_1rhs
        procedure:: solveLinearLeastSquare_nrhs
        procedure:: solveLinearLeastSquare_Constrain
    end interface solveLinearLeastSquare
    
contains

    
    !solve the symmetry linear equations system by LU factorization
    pure subroutine solveSymmetryLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::   a
    real(rp),dimension(1:),intent(inout)::      b
    integer(isp),allocatable,dimension(:)::     ipiv
        !refer to ?getrf of lapack
        allocate( ipiv( max(1,min(size(a,dim=1),size(a,dim=2))) ) )
        call sytrf(a,ipiv=ipiv)
        call sytrs(a,b,ipiv)
    end subroutine solveSymmetryLES
    
    
    ![d] is the diagonal vector with size(n)
    ![eva] input sub_diagonal in beginning n-1 position and output n eigenvalues
    ![evc] output n*n orthonomal eigenvector
    !give the eigenvalues and eigenvectors for symmetric Tridiagnoal matrix
    pure subroutine eigenSymTriDiagonal(d,eva,evc)
    real(rp),dimension(:),intent(in)::              d
    real(rp),dimension(:),intent(inout)::           eva
    real(rp),dimension(:,:),intent(out),optional::  evc
        if(present(evc)) then
            call stev(d,eva,evc)
        else
            call stev(d,eva)
        endif
    end subroutine eigenSymTriDiagonal
    
    
    !call steqr(d,e,z,compz,info)
    ![d]    input the diagonal vector with size(n) and output the n eigenvalues in ascending order
    ![e]    input contains the off-diagonal elements of T with size(n-1)
    ![evc]  output n*n orthonomal eigenvector
    !Find all eigenvalues and eigenvectors of a tridiagonal matrix T
    pure subroutine eigenTriDiagonal(d,e,evc)
    real(rp),dimension(:),intent(inout)::               d
    real(rp),dimension(:),intent(in)::                  e
    complex(rp),dimension(:,:),intent(out),optional::   evc
        if(present(evc)) then
            call steqr(d,e,evc)
        else
            call steqr(d,e)
        endif
    end subroutine eigenTriDiagonal
    
    
    !dir$ if defined (noMKL)
    !because triFacotrSquareMat may lead U be singluar, this procedure may fails
    pure subroutine solveGeneralLES(a,b)
    real(rp),dimension(:,:),intent(inout)::     a
    real(rp),dimension(:),intent(inout)::       b
    real(rp),dimension(size(b))::               y
    integer(ip)::                               n,i
        call triFactorSquareMat(a)
        n = size(b)
        !solve LY=b
        y(1) = b(1)
        do i=2,n
            y(i) = b(i) - sum( a(i,1:i-1) * y(1:i-1) )
        enddo
        !solve UX=Y where b as X
        b(n) = y(n) / a(n,n)
        do i=n-1,1,-1
            b(i) = ( y(i) - sum( a(i,i+1:n) * b(i+1:n) ) ) / a(i,i)
        enddo
    end subroutine solveGeneralLES
    !dir$ else
    !refer to ?getrf of lapack
    !getrf is short for triangular factorization of general real matrices
    !solve the general linear equations system by triangular factorization
    pure subroutine solveGeneralLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::   a
    real(rp),dimension(1:),intent(inout)::      b
    integer(isp),allocatable,dimension(:)::     ipiv
        allocate(ipiv(max(1 , min(size(a,dim=1) , size(a,dim=2)))))
        call getrf(a,ipiv)
        call getrs(a,ipiv,b)
    end subroutine solveGeneralLES
    !dir$ end if
    
    
    !dir$ if .not. defined(noMKL)
    !solve general linear least square
    !m number of rows | n numer of columns | nrhs number of right-hand side, number of columns in B
    !input: a(m,n) | b(max(m,n) , nrhs)
    !output: a QR factorization if(m>=n) and LQ factorization if(m<n)
    !output: b solution vector (1:n,nrhs)
    !m > n overdetermined | solve for least square
    !m < n underdetermined | solve for minimum norm solution
    pure subroutine solveLinearLeastSquare_1rhs(a,b)
    real(rp),dimension(:,:),intent(inout):: a
    real(rp),dimension(:),intent(inout)::   b
    real(rp),dimension(size(b),1)::         b0
        b0(:,1) = b
        call gels(a,b0)
        b = b0(:,1)
    end subroutine solveLinearLeastSquare_1rhs

    !see https://www.quora.com/Is-it-better-to-do-QR-Cholesky-or-SVD-for-solving-least-squares-estimate-and-why
    !the different least square methods based on QR or SVD 
    pure subroutine solveLinearLeastSquare_nrhs(a,b)
    real(rp),dimension(:,:),intent(inout):: a,b
        call gels(a,b)  !QR
        !call gelss(a,b) !SVD, sigular value decomposition
    end subroutine solveLinearLeastSquare_nrhs

    !solve min|c-ax| s.t. bx=d
    !the dimension: a(m,n) | c(m) | b(p,n) | d(p) | x(n)
    !should be: p<=n<=m+p; rank(b)=p; rank(a/b)=n (full rank)
    pure subroutine solveLinearLeastSquare_Constrain(a,b,c,d,x,info)
    real(rp),dimension(:,:),intent(inout):: a,b
    real(rp),dimension(:),intent(inout)::   c,d
    real(rp),dimension(:),intent(out)::     x
    integer(ip),optional,intent(out)::      info
    integer(4)::                            i
        call gglse(a,b,c,d,x,i)
        if(present(info)) info = i
    end subroutine solveLinearLeastSquare_Constrain
    !dir$ endif
    
    
    !dir$ if defined (noMKL)
    !dir$ else
    pure subroutine inverseGeneralSquareMat(a)
    real(rp),dimension(:,:),intent(inout)::     a
    integer(isp),allocatable,dimension(:)::     ipiv
        allocate(ipiv(max(1 , min(size(a,dim=1) , size(a,dim=2)))))
        call getrf(a,ipiv)
        call getri(a,ipiv)
    end subroutine inverseGeneralSquareMat
    !dir$ end if
    
    !dir$ if defined (noMKL)
    !dir$ else
    !if a is positive-definite matrix, it can be factor as a=R^TR
    !in this subroutine, we assume a is symmetric
    pure subroutine choleskyFactorSym(a)
    real(rp),dimension(:,:),intent(inout)::     a
        call potrf(a,'U')
    end subroutine choleskyFactorSym
    !dir$ end if
    
    !dir$ if defined (noMKL)
    !input a [(2:n,1),(1:n,2),(1:n-1,3)]
    !refer to chasing method
    !limiting: abs(a)=>(a(1,2)>a(1,3)>0),((a(i,2)>a(i,1)+a(i,3)),(a(n,2)>a(n,1))
    !a validation refer to 
    !https://wenku.baidu.com/view/a2065cb064ce0508763231126edb6f1aff0071d7.html
    !the lapackwrapper and chasing method seem no difference below...
    pure subroutine solveTridiagonalLES(a,b)
    real(rp),dimension(1:,1:),intent(in)::  a
    real(rp),dimension(1:),intent(inout)::  b
    real(rp),dimension(:),allocatable::     beta
    integer(ip)::                           i,n
        n = size(a,dim=1)
        allocate(beta(n-1))
        !step1
        beta(1) = a(1,3) / a(1,2)
        do i=2,n-1
            beta(i) = a(i,3) / ( a(i,2) - a(i,1) * beta(i-1) )
        enddo
        !step2
        b(1) = b(1) / a(1,2)
        do i=2,n
            b(i) = ( b(i) - a(i,1) * b(i-1) ) / ( a(i,2) - a(i,1) * beta(i-1) )
        enddo
        !step3
        b(n) = b(n)
        do i=n-1,1,-1
            b(i) = b(i) - beta(i) * b(i+1)
        enddo
    end subroutine solveTridiagonalLES
    !dir$ else    
    pure subroutine solveTridiagonalLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::a
    real(rp),dimension(1:),intent(inout)::  b
    integer(ip)::                           n
        n = size(a,dim=1)
        call dttrfb(a(2:n,1),a(1:n,2),a(1:n-1,3))
        call dttrsb(a(2:n,1),a(1:n,2),a(1:n-1,3),b)
    end subroutine solveTridiagonalLES
    !dir$ end if
    
    

!--------------------------------------------------------------
    !emergency only and very limited
    !m[in] is matrix, and m[out] is [L-I+U] witout the unit dialog of L
    !this is only for square matrix, and without ipiv unlike getrf in lapack
    !--
    !this simple procedure may fail, details can refer to <LU decompostion> in wiki
    !considering a11=0, then u11=0, U is singular. but A may be not singular, so a logic error 
    !a Partial Pivoting is needed for a proper implementation, and waiting for improving
    pure subroutine triFactorSquareMat(m)
    real(rp),dimension(:,:),intent(inout):: m
    integer(ip)::                           i,j,n
    
        !dir$ if defined (lwcheck)
        if(size(m,dim=1)/=size(m,dim=2)) call disableProgram
        !dir$ end if

        n = size(m,dim=1)
        m(2:n,1) = m(2:n,1) / m(1,1)
        do j=2,n-1
            !upper part | if(i==1) m(i,j) = m(i,j)
            do i=2,j
                m(i,j) = m(i,j) - sum(m(i,1:i-1)*m(1:i-1,j))
            enddo
            !lower part
            do i=j+1,n
                m(i,j) = ( m(i,j) - sum(m(i,1:j-1)*m(1:j-1,j)) ) / m(j,j)
            enddo
        enddo
        do i=2,n
            m(i,n) = m(i,n) - sum(m(i,1:i-1)*m(1:i-1,n))
        enddo
    
    end subroutine triFactorSquareMat
    
    
    
!-------------------------------------------------
    !gmres: suppose x is in the space x_0 + \mathcal(K)_m, then we have expression
    !x = x_0 + V_m y
    !to solve x for min{J = |b - A x|}, then transfer to
    !x_m = x0 + V_m y_m
    !y_m = argmin_y | \beta e_1 - \hat{H}_m y |_2
!-------------------------------------------------
    !solve AX=b, A is sparse matrix, full size (n*n)
    ![nnz] means number of non-zero
    !csr structure means compressed sparse row, consists of three parts(a(nnz),ia(n+1),ja(nnz))
    !a: the value of non-zero
    !ia: the start of a row of the sparse matrix in <nnz> dimension, range(1:nnz)
    !ja: the location in this row, range (1:n)
    !maxInr range (1:n)
    subroutine gmres(a,ia,ja,b,maxOtr,maxInr, x)
    real(rp),dimension(:),intent(in)::      a,b
    integer(ip),dimension(:),intent(in)::   ia,ja
    integer(ip),intent(in)::                maxOtr,maxInr
    real(rp),dimension(:),intent(inout)::   x
    real(rp),parameter::                    rEps = 10._rp * GlobalEps
    real(rp),dimension(size(b))::           r
    real(rp),dimension(maxInr)::            c,s
    real(rp),dimension(maxInr+1)::          g,y
    real(rp),dimension(size(b),maxInr+1)::  v
    real(rp),dimension(maxInr+1,maxInr)::   h
    integer(ip)::                           i,j,iter,jcopy
    real(rp)::                              beta,eps,av,t,mu
    
        !outer loop
        do iter = 1,maxOtr
            
            call csrAx(a,ia,ja,x,r); r = b - r
            beta = norm2(r)
            if(iter==1) eps = beta * rEps
            v(:,1) = r / beta
            
            g = 0._rp; g(1) = beta
            h = 0._rp
            
            !inner loop
            do j=1,maxInr
            
                jcopy = j
            
                call csrAx(a,ia,ja,v(:,j),v(:,j+1))
                av = norm2(v(:,j+1))
                
                !--orth
                do i=1,j
                    h(i,j) = v(:,j+1) .ip. v(:,i)
                    v(:,j+1) = v(:,j+1) - h(i,j)*v(:,i)
                enddo
                h(j+1,j) = norm2(v(:,j+1))
                
                !special deal for singular sparse matrix
                if(h(j+1,j)/av < rEps) then
                    do i=1,j
                        t = v(:,j+1) .ip. v(:,i)
                        h(i,j) = h(i,j) + t
                        v(:,j+1) = v(:,j+1) - t*v(:,i)
                    enddo
                    h(j+1,j) = norm2(v(:,j+1))
                endif
                
                if(h(j+1,j)/=0._rp) v(:,j+1) = v(:,j+1) / h(j+1,j)
                
                !rot
                if(j>1) then
                    do i=1,j-1
                        h(i:i+1,j) = rot2([c(i),s(i)],h(i:i+1,j))
                    enddo
                endif
                
                mu = sqrt(sum(h(j:j+1,j)**2))
                c(j) = h(j,j)/mu
                s(j) = -h(j+1,j)/mu
                h(j,j) = c(j) * h(j,j) - s(j) * h(j+1,j)
                h(j+1,j) = 0._rp
                
                g(j:j+1) = rot2([c(j),s(j)] , g(j:j+1))
                beta = abs(g(j+1))
                
                if(beta < eps) exit
                
            enddo
            
            j = jcopy - 1
            y(j+1) = g(j+1)/h(j+1,j+1)
            
            do i=j,1,-1
                y(i) = (g(i) - (h(i,i+1:j+1) .ip. y(i+1:j+1))) / h(i,i)
            enddo
            
            do i=1,size(b)
                x(i) = x(i) +  (v(i,1:j+1) .ip. y(1:j+1))
            enddo
            
            if(beta < eps) exit
            
        enddo
        
    contains
    
        !--
        pure subroutine Arnoldi(v,H,ksp)
        real(rp),dimension(:),intent(in)::      v
        real(rp),dimension(:,:),intent(inout):: H,ksp
        
            
        
        end subroutine Arnoldi
        
        
        !--
        pure subroutine HouseholderArnoldi()
        !real(rp),dimension(:),intent(in)::      
        
        end subroutine HouseholderArnoldi
        
    end subroutine gmres

!-----------------------------------------------    
    !dir$ if .not. defined(noMKL)
    !Ax=b | A=A(a,ia,ja) | input x0 output x |
    subroutine fgmresILUT(a,ia,ja,b,x)
    real(8),dimension(:),intent(in)::       a
    integer(4),dimension(:),intent(in)::    ia,ja
    real(8),dimension(:),intent(inout)::    x,b
    real(8),dimension(size(x))::            trv
    integer(4)::                            n,nm1,rciRequest,iter,ip15,maxfil,ierr,itercount
    integer(4),dimension(128)::             ipar
    integer(4),dimension(:),allocatable::   ibilut,jbilut
    real(8)::                               tol
    real(8),dimension(128)::                dpar
    real(8),dimension(:),allocatable::      tmp,bilut
    
        n = size(x); nm1 = n-1
        !--
        maxfil = 1  !half bandwidth of preconditioner M=LU
        allocate(bilut((2*maxfil+1)*n - maxfil*(maxfil+1) + 1))
        allocate(jbilut((2*maxfil+1)*n - maxfil*(maxfil+1) + 1))
        allocate(ibilut(n+1))
        
        !--
        ip15 = min(2,n) !restarted gmres, default min(150,n) which means no restarted gmres
                        !sometimes small ip15 is sufficient for convergence and sometimes large ipar15 is necessary
        allocate(tmp((2*ip15+1)*n + pidb2(ip15*(ip15+9)) + 1))
        
        !--
        call dfgmres_init(n, x, b, rciRequest, ipar, dpar, tmp)
        
        call preconditioner
        
        call setCheck
        
        !loop
        do; call dfgmres(n,x,b,rciRequest,ipar,dpar,tmp)
            select case(rciRequest)
            case(0)
                exit
            case(1)
                call csrAx(a,ia,ja, tmp(ipar(22):ipar(22)+nm1), tmp(ipar(23):ipar(23)+nm1))
            case(2) !stopping test
                stop 'error: laWrapperLib/fgmresILUT has no user-defined stopping test, ipar(10)=0'
            case(3)
                call mkl_dcsrtrsv('l','n','u',n,bilut,ibilut,jbilut,tmp(ipar(22)),trv)
                call mkl_dcsrtrsv('u','n','n',n,bilut,ibilut,jbilut,trv,tmp(ipar(23)))
            case(4)
                if(dpar(7)<1.e-12_8) exit   !norm of generated vector
            case default
                print*, 'fgmres request:  ', rciRequest
                stop 'error: laWrapperLib/fgmresILUT fails'
            end select
        enddo
        
        ipar(13) = 0 !update x
        call dfgmres_get(n, x, b, rciRequest, ipar, dpar, tmp, itercount)
        call mkl_free_buffers
        
    contains
        !--
        subroutine preconditioner
            ipar(31) = 1        !0:abort ilu of routine meets zero diag element; 1: plus dpar(31)*norm(ev)
            dpar(31) = 1.e-7_8  !diag value
            tol = 1.e-8_8       !tolerance
            call dcsrilut(n, a, ia, ja, bilut, ibilut, jbilut, tol, maxfil, ipar, dpar, ierr)
        end subroutine preconditioner
        !--
        subroutine SetCheck
            ipar(8) = 0         !not perform stop testing for loop | ipar(4)[current outer loop] < ipar(5)[min(150,n)]
            ipar(9) = 1         !not perform stop testing for accuracy | dpar(5)[currentResNorm] < dpar(4)
                                != dpar(1)*dpar(3)_0 + dpar(2) [relativeTol * initResNorm + absoluteTol]
            ipar(10) = 0        !do perform user-defined stopping test
            ipar(11) = 1        !0: run non-preconditioned dfgmres; 1: run preconditioned dfgmres
            ipar(12) = 0        !not perform norm of currently generated vecotr | dpar(7) < dpar(8) [1.e-12_8]
            ipar(15) = ip15     !restarted gmres[inner loop], default min(150,n) which means no restarted gmres
            dpar(1) = 1.e-8_8
            call dfgmres_check(n, x, b, rciRequest, ipar, dpar, tmp)
        end subroutine setCheck
    end subroutine fgmresILUT
    
    subroutine fgmres(a,ia,ja,b,x)
    real(8),dimension(:),intent(in)::       a
    integer(4),dimension(:),intent(in)::    ia,ja
    real(8),dimension(:),intent(inout)::    x,b
    integer(4)::                            n,nm1,rciRequest,ip15,itercount
    integer(4),dimension(128)::             ipar
    real(8),dimension(128)::                dpar
    real(8),dimension(:),allocatable::      tmp
    
        n = size(x); nm1 = n-1
        !--
        ip15 = min(2,n) !restarted gmres, default min(150,n) which means no restarted gmres
                        !sometimes small ip15 is sufficient for convergence and sometimes large ipar15 is necessary
        allocate(tmp((2*ip15+1)*n + pidb2(ip15*(ip15+9)) + 1))

        !--
        call dfgmres_init(n, x, b, rciRequest, ipar, dpar, tmp)
        
        call setCheck
        
        !loop
        do; call dfgmres(n,x,b,rciRequest,ipar,dpar,tmp)
            select case(rciRequest)
            case(-1)
                print*, 'warning: laWrapperLib/fgmres get a bad initial guess, solution may be inaccurate'
                exit !
            case(0)
                exit
            case(1)
                call csrAx(a,ia,ja, tmp(ipar(22):ipar(22)+nm1), tmp(ipar(23):ipar(23)+nm1))
            case(2) !stopping test
                stop 'error: laWrapperLib/fgmres has no user-defined stopping test, ipar(10)=0'
            case(3)
                stop 'error: laWrapperLib/fgmres has no user-defined preconditioner'
            case(4)
                stop 'error: laWrapperLib/fgmres has no user-defined vector check'
            case default
                print*, 'fgmres request:  ',rciRequest
                stop 'error: laWrapperLib/fgmres fails'
            end select
        enddo
        
        ipar(13) = 0 !update x
        call dfgmres_get(n, x, b, rciRequest, ipar, dpar, tmp, itercount)
        call mkl_free_buffers
        
    contains
        !--
        subroutine SetCheck
            ipar(8) = 1         !perform stop testing for loop | ipar(4)[current outer loop] < ipar(5)[min(150,n)]
            ipar(9) = 1         !perform stop testing for accuracy | dpar(5)[currentResNorm] < dpar(4)
                                != dpar(1)*dpar(3) + dpar(2) [relativeTol * initResNorm + absoluteTol]
            ipar(10) = 0        !do not perform user-defined stopping test
            ipar(11) = 0        !0: run non-preconditioned dfgmres; 1: run preconditioned dfgmres
            ipar(12) = 1        !perform norm of currently generated vecotr | dpar(7) < dpar(8) [1.e-12_8]
            ipar(15) = ip15     !restarted gmres[inner loop], default min(150,n) which means no restarted gmres
            dpar(1) = 1.e-8_8   !relative tolerance
            dpar(2) = 0._8      !abolute tolerance
            call dfgmres_check(n, x, b, rciRequest, ipar, dpar, tmp)
        end subroutine setCheck
    end subroutine fgmres
    !dir$ else
    
    !dir$ end if
    
    
    !dir$ if .not. defined(noMKL)
    !ge[sparse representation of a general matrix]
    !mv[matrix-vector product ]
    subroutine csrAx(a,ia,ja,x,b)
    real(rp),dimension(:),intent(in)::      a,x
    integer(ip),dimension(:),intent(in)::   ia,ja
    real(rp),dimension(:),intent(out)::     b
        call mkl_dcsrgemv('n', size(x), a, ia, ja, x, b)
    end subroutine csrAx
    !dir$ else
    pure subroutine csrAx(a,ia,ja,x,b)
    real(rp),dimension(:),intent(in)::      a,x
    integer(ip),dimension(:),intent(in)::   ia,ja
    real(rp),dimension(:),intent(out)::     b
    integer(ip)::                           i,j,k
        b = 0._rp
        do i=1,size(x)
            j = ia(i);
            k = ia(i+1) - 1
            b(i) = b(i) + (a(j:k) .ip. x(ja(j:k)))
        enddo
    end subroutine csrAx
    !dir$ end if
    
    pure function polyfit(x,y,n)
    real(rp),dimension(:),intent(in)::  x,y
    integer(ip),intent(in)::            n
    real(rp),dimension(0:n)::           polyfit
    real(rp),dimension(size(x),n+1)::   A
    real(rp),dimension(size(x))::       yt
    integer(ip)::                       m,i,j
    
        m = size(x)
        if(m <= n) call disableprogram !underdetermination

        yt = y
        do j=1,n+1
            do i=1,m
                A(i,j) = x(i)**(j-1)
            enddo
        enddo
        
        call solveLinearLeastSquare(A, yt)
        polyfit = yt(1:n+1)
        
    end function polyfit
    
end module laWrapperLib