!this is a space(mathematics) based on the askey polynomial basis
!a function space with different measures has different optimal basis
!it is a space due to the commutative operation (+ and *) and the measurement induced by inner product
!the operational object is the coordinates of truncated basis expressed as an array
module polynomialSpace_
use constants
use arrayOpsLib
use stringOpsLib
use IntegrationLib
use laWrapperLib
use quadprogWrapperLib
use polynomial_
use specialPolynomial_
implicit none

    private
    public:: polynomialSpace
    public:: polyType
    
    !--
    character(10),dimension(6),parameter:: polyType = ['legendre','hermite',    &
            'jacobi','chebyshev1','chebyshev2','laguerre']
    
!----------------------------------
    !due to the normal-orthogonal polynomial basis, assert: basis(0) = 1 => phi(0) = mean(phi)
    !details see my TexNotes/UQ/fragment/method
    type:: polynomialSpace
    
        private
        
        integer(ip)::                               quadnp_
        integer(ip)::                               truncOd_
        
        !--
        type(polynomial),dimension(:),allocatable:: basis_
        character(10)::                             basisType_
        !the innerproduct of the space is based on the quadrature rule
        !innerproduct and quadarture rule may differ with a coef due to 
        !the normalized measurement/probability
        real(rp)::                                  ipMeasCoef_
        
        !(0:so,1:np), use quadrature rule for innerproduct
        real(rp),dimension(:,:),allocatable::       quadBasisVal_
        real(rp),dimension(:),allocatable::         quadx_
        real(rp),dimension(:),allocatable::         quadw_
        
        !--A a = funcval(x_gp) | A = trans(quadBasisVal_)
        real(rp),dimension(:,:),allocatable::       ata_
        
        !use tribasis quadrature value for multiplication and division
        !triMat_{ijk} = \int \phi_i \phi_j \phi_k d(P(\xi))
        real(rp),dimension(:,:,:),allocatable::     triMat_
        
    contains
    
        generic::               init => init_ba, init_od
        procedure::             init_ba
        procedure::             init_od
        
        procedure::             makeOpsCoef

        !--
        generic::               quadrule => quadrule_unit, quadrule_multi
        procedure::             quadrule_unit
        procedure::             quadrule_multi
        
        !--
        generic::               quadval => quadval_x, quadval_jox, quadval_ig, quadval_igjo
        procedure::             quadval_x       !return function val at x
        procedure::             quadval_jox     !return jth order basis val at point x
        procedure::             quadVal_ig      !return function val at gauss point i
        procedure::             quadval_igjo    !return jth order basis val at gauss point i
        
        !--
        generic::               project => project_poly, project_func1
        procedure::             project_poly    !function of polynomial form 
        procedure::             project_func1   !one dimensional function
        
        !--
        generic::               limit => limit_linear_positive
        procedure::             limit_linear_positive
        
        generic::               limit2 => limit_quadprog_positive
        procedure::             limit_quadprog_positive
        
        !--
        procedure::             writefunc => writefunc_ps
        
        !operation in space, short for use
        !--binary
        procedure::             mt                  !multiply
        procedure::             dv                  !divide
    
        !--unitary
        generic::               op => op1,op2
        procedure::             op1                 !user-defined unitary operator
        procedure::             op2                 !user-defined binary operator
    
        procedure::             sqt                 !sqrt
        procedure::             ex                  !exponent
        procedure::             ln                  !log_e
        generic::               pw => pw_i, pw_r    !power
        procedure::             pw_i, pw_r
        
        generic::               dc => dc_func, dc_val   !discontinuity
        procedure::             dc_func
        procedure::             dc_val
        
        !--derivative
        !\frac{\partial f_j}{\partial x_k} = \int \frac{\partial f}{\partial x} \frac{\partial x}{\partial x_k} 
        !\phi_j = \sum (\frac{\partial f}{\partial x})_i S_ijk
        procedure::             deri    !derivative
        !diff = df = df/du du 
        procedure::             diff    !differential
        
        !member function
        procedure::             truncOd !
        procedure::             quadnp  !
        procedure::             quadx   !
        procedure::             triMat
    
    end type polynomialSpace
    
contains

    subroutine init_ba(this, basisType, basis, np)
    class(polynomialSpace),intent(out)::        this
    character(*),intent(in)::                   basisType
    type(polynomial),dimension(0:),intent(in):: basis
    integer(ip),optional,intent(in)::           np
    integer(ip)::                               truncOrder,i
    character(len(basisType))::                 T
    
        truncOrder = ubound(basis,1)
        this%truncOd_ = truncOrder
        
        !--3*so/2+1 is safe for triBasis, but not enough for other operator
        !quadnp>0 when truncorder=0
        this%quadnp_ = merge(np, 4*(truncOrder+1), present(np))
        
        !--
        allocate(this%basis_(0:ubound(basis,1)), source = basis)
        
        !--
        t = basisType
        call lowerString(t)
        call check
        this%basisType_ = t
        
        !--
        select case(t)
        case(polyType(1))
            this%ipMeasCoef_ = 1._rp/2._rp
        case(polyType(2))
            this%ipMeasCoef_ = 1._rp/sqrt(2._rp*pi)
        end select
        
        !--a basisType alwasy correspond to a quadrature type
        call this%makeOpsCoef(basis, basisType, this%quadnp_, this%ipMeasCoef_)
        
    contains
    
        subroutine check
        
            do i=1,size(polyType)
                if(t == polyType(i)) exit
            end do
            
            if(i==size(polyType)+1) &
            stop 'error: polynomialSpace get an unrecognized type'
            
            if(truncOrder<0) &
            stop 'error: polynomialSpace get a negative order'
            
        end subroutine check
        
    end subroutine init_ba
    
    !--
    subroutine init_od(this,basisType,od,np)
    class(polynomialSpace),intent(out)::    this
    character(*),intent(in)::               basisType
    integer(ip),intent(in)::                od
    integer(ip),optional,intent(in)::       np
    character(len(basisType))::             t
    type(polynomial),dimension(0:od)::      ba
    
        t = basisType
        call lowerString(t)
        
        select case(t)
        case(polytype(1))
            ba = normalLegendrepolynomialSet(od)
            ba = sqrt(2._rp)*ba
        case(polytype(2))
            ba = normalHermitePolynomialSet(od, 'prob')
            ba = sqrt(sqrt(2._rp)*srpi)*ba
        case default
            stop 'error: polynomialSpace/init_od get error type'
        end select

        if(present(np)) then
            call this%init(t, ba, np)
        else
            call this%init(t, ba)
        endif
    
    end subroutine init_od
    
    
    !--
    subroutine makeOpsCoef(this,basis,quadtype,quadnp,ipMeasCoef)
    class(polynomialSpace),intent(inout)::      this
    type(polynomial),dimension(0:),intent(in):: basis
    character(*),intent(in)::                   quadtype
    integer(ip),intent(in)::                    quadnp
    real(rp),intent(in)::                       ipMeasCoef
    real(rp)::                                  x
    integer(ip)::                               od,np,ncnp,i,j,k
        
        if(allocated(this%quadw_)) deallocate(this%quadw_)
        if(allocated(this%quadBasisVal_)) deallocate(this%quadBasisVal_)
        if(allocated(this%triMat_)) deallocate(this%triMat_)
        np = quadnp
        od = ubound(basis, 1)
        
        !assume the accuracy of quadrature rule is N(clenshawcurtis eg.), 
        !guarantee the tribasis computation
        if(np<3*od) stop 'error: AskyPolynomialSpace/makeOpsCoef get a bad np or order'
        
        !do not store x, but the value of polynomial in x
        allocate(this%quadx_(np) , this%quadw_(np))
        if(quadtype==polytype(2)) then
            call quadratureRule('hermite_prob', this%quadx_, this%quadw_)
        else
            call quadratureRule(quadtype, this%quadx_, this%quadw_)
        endif

        !---
        allocate(this%quadBasisVal_(0:od,1:np))
        do j=1,np
            do i=0,od
                this%quadBasisVal_(i,j) = basis(i)%funcval(this%quadx_(j))
            enddo
        enddo
        
        !--A = trans(quadBasisVal_)
        allocate(this%ata_, source = this%quadBasisVal_.ip.trans(this%quadBasisVal_))
        
        !tribasis method is more robust for multiplication and division than quadrature method
        allocate(this%triMat_(0:od, 0:od, 0:od))
        do k=0,od
            do j=0,k
                do i=0,j
                    this%triMat_(i,j,k) = ipMeasCoef*sum( &
                        (this%quadBasisVal_(i,:)*this%quadBasisVal_(j,:)*this%quadBasisVal_(k,:)) &
                        *this%quadw_)
                    this%triMat_(i,k,j) = this%triMat_(i,j,k)
                    this%triMat_(j,i,k) = this%triMat_(i,j,k)
                    this%triMat_(j,k,i) = this%triMat_(i,j,k)
                    this%triMat_(k,i,j) = this%triMat_(i,j,k)
                    this%triMat_(k,j,i) = this%triMat_(i,j,k)
                enddo
            enddo
        enddo
        
    end subroutine makeOpsCoef
    
    
    !--ipmc: measure coefficients for inner product
    subroutine quadrule_unit(this, x, w)
    class(polynomialspace),intent(in)::     this
    real(rp),dimension(:),intent(out)::     x,w
        if(this%basisType_==polytype(2)) then
            call quadratureRule('hermite_prob', x, w)
        else
            call quadratureRule(this%basisType_, x, w)
        endif
    end subroutine quadrule_unit
    
    !--
    subroutine quadrule_multi(this, x, w)
    class(polynomialspace),intent(in)::     this
    real(rp),dimension(:,:),intent(out)::   x
    real(rp),dimension(:),intent(out)::     w
        x = 0._rp
        w = 0._rp
        stop 'error: polynomialSpace for quadrule_multi is bulding'
    end subroutine quadrule_multi
    

    !============================================================
    elemental real(rp) function quadval_igjo(this,quadi,orderj)
    class(polynomialSpace),intent(in):: this
    integer(ip),intent(in)::            quadi,orderj
        quadval_igjo = this%quadBasisVal_(orderj, quadi)
    end function quadval_igjo
    
    !return the value at the quad point based on the coefficients[a]
    pure real(rp) function quadval_ig(this,a,quadi)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    integer(ip),intent(in)::            quadi
        quadval_ig = sum(a*this%quadBasisVal_(:, quadi))
    end function quadval_ig
    
    pure real(rp) function quadval_x(this,a,x)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),intent(in)::               x
        quadval_x = sum(a*this%basis_([0:this%truncOd()])%funcval(x))
    end function quadval_x
    
    elemental real(rp) function quadval_jox(this,orderj,x)
    class(polynomialSpace),intent(in):: this
    integer(ip),intent(in)::            orderj
    real(rp),intent(in)::               x
        quadval_jox = this%basis_(orderj)%funcval(x)
    end function quadval_jox
    
    !====
    pure real(rp) function project_func1(this,f1,ibasis) result(c)
    class(polynomialSpace),intent(in):: this
    procedure(absf1)::                  f1
    integer(ip),intent(in)::            ibasis
        c = this%ipMeasCoef_*sum(f1(this%quadx_)*this%quadBasisVal_(ibasis,:)*this%quadw_)
    end function project_func1
    
    elemental real(rp) function project_poly(this,polyfunc,ibasis) result(c)
    class(polynomialSpace),intent(in):: this
    type(polynomial),intent(in)::       polyfunc
    integer(ip),intent(in)::            ibasis
        c = this%project(pv, ibasis)
    contains
        elemental real(rp) function pv(x)
        real(rp),intent(in)::   x
            pv = polyfunc%funcval(x)
        end function pv
    end function project_poly
    
    !-----
    subroutine limit_linear_positive(this,a)
    class(polynomialspace),intent(in)::     this
    real(rp),dimension(0:),intent(inout)::  a
    integer(ip)::                           i
    real(rp)::                              quadval,beta,b,eps
    
        beta = minrp
        eps = GlobalEps*a(0)
        do i=1,this%quadnp_
            !the val at quadrature point based on the spectrum a
            quadval = this%quadval(a, i)
            if(quadval<0._rp) then
                !this b make u0 + (1-b)*sum(u_i\phi_i) = eps
                b = quadval/(quadval - a(0)) + eps
                b = max(0._rp, min(b, 1._rp))
                !leave the max b which guaratee the function at all of quadrature point positive
                if(b>beta) beta = b
            endif
        enddo
        
        if(beta>0._rp) a(1:this%truncOd_) = (1-beta)*a(1:this%truncOd_)
        
    end subroutine limit_linear_positive
    
    !limit_quadprog_positive
    subroutine limit_quadprog_positive(this,a)
    class(polynomialspace),intent(in)::     this
    real(rp),dimension(0:),intent(inout)::  a
    integer(ip)::                           i, n, q, meq
    real(rp),dimension(:,:),allocatable::   dmat, amat
    real(rp),dimension(:),allocatable::     d,b
    
        do i=1,this%quadnp_
            if(this%quadval(a, i)<0._rp) exit
        enddo
        if(i==this%quadnp_+1) return
        
        !n is the dimension of quadrature equations
        !q is the dimension of limitation    
        n = size(a)
        q = this%quadnp_
        allocate(dmat(n,n), d(n), amat(n,q), b(q))
        
        !min( |A[q,n]a[n] - b[q]|_2 ) | quandbasisval_(0:so,1:np) | n=so+1, q=np
        !A = trans(this%quadBasisVal_)
        !D = A^T A | this%quadBasisVal_ .ip. trans(this%quadBasisVal_)
        dmat(:,:) = this%ata_
        !d^T = b^T A = (A a)^T A = a A^T A
        d(:) = a .ip. this%ata_
        !--
        amat(:,:) = this%quadBasisVal_
        b(:) = GlobalEps*a(0)
        
        meq = 0
        call quadprog(dmat, d, amat, b, meq, a)
        
    end subroutine limit_quadprog_positive
    
    !------------------------
    subroutine writefunc_ps(this,a,lo,up,np,filename)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),intent(in)::               lo,up
    integer(ip),intent(in)::            np
    character(*),optional,intent(in)::  filename
    
        if(present(filename)) then
            call writefunc(val, lo, up, np, filename)
        else
            call writefunc(val, lo, up, np, 'psfunc')
        endif
        
    contains
    
        elemental real(rp) function val(xi)
        real(rp),intent(in)::   xi
            val = this%quadval(a,xi)
        end function val
    
    end subroutine writefunc_ps
    
    
    !==========================================
    !operator
    !==========================================
    pure function mt(this,lhs,rhs)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: lhs,rhs
    real(rp),dimension(0:ubound(lhs,1)):: mt
    integer(ip)::                       i,j,k
        mt = 0._rp
        do k=0,this%truncOd_
            do j=0,this%truncOd_
                do i=0,this%truncOd_
                    mt(k) = mt(k) + lhs(i)*rhs(j)*this%triMat_(i,j,k)
                enddo
            end do
        enddo
    end function mt
    
    !--
    pure function dv(this,up,down)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: up,down
    real(rp),dimension(0:ubound(down,1)):: dv
    real(rp),dimension(0:ubound(down,1),0:ubound(down,1)):: mat
    integer(ip)::                       i,j
        do j=0,this%truncOd_
            do i=0,this%truncOd_
                mat(i,j) = sum(down(:)*this%triMat_([0:this%truncOd_], j, i))
            enddo
        enddo
        dv = up
        call solveGeneralLES(mat, dv)
    end function dv
    
    !--unitary operation o_k = \int func(u) \phi_k dp(\xi) = \sum func(\xi_i) \phi_k(\xi_i) w(i)
    !--if let func = df/du = df/du(u), it can be used to construct the derivative
    !--df^i/du_j = this%deri(this%op(a,df/du)) where u = \sum a_i \phi_i
    pure function op1(this,a,func) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    procedure(absf1)::                  func
    real(rp),dimension(0:ubound(a,1)):: o
    real(rp),dimension(:),allocatable:: quadKernel
    integer(ip)::                       i,n,np
    
        n = this%truncOd_
        np = this%quadnp_
        allocate(quadKernel(np))
        do i=1,np
            quadKernel(i) = func(this%quadval(a, i))*this%quadw_(i)
        enddo
        do i=0,n
            o(i) = this%ipMeasCoef_*sum(quadKernel*this%quadBasisVal_(i,:))
        enddo
        
    end function op1
    
    pure function op2(this,a,b,func) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a,b
    procedure(absf2)::                  func
    real(rp),dimension(0:ubound(a,1)):: o
    real(rp),dimension(:),allocatable:: quadKernel
    integer(ip)::                       i,n,np
    real(rp)::                          x,y
        n = this%truncOd_
        np = this%quadnp_
        allocate(quadKernel(np))
        do i=1,np
            x = this%quadval(a, i)
            y = this%quadval(b, i)
            quadKernel(i) = func(x,y)*this%quadw_(i)
        enddo
        do i=0,n
            o(i) = this%ipMeasCoef_*sum(quadKernel*this%quadBasisVal_(i,:))
        enddo
    end function op2
    
    !--
    pure function sqt(this,a) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),dimension(0:ubound(a,1)):: o
        !solve \int \sqrt(\sum_{k=0}^{n} u_k*\psi_k) \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = sqrt(x)
        end function ker
    end function sqt
    
    !--
    pure function pw_r(this,a,s) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),intent(in)::               s
    real(rp),dimension(0:ubound(a,1)):: o
        !solve \int (\sum_{k=0}^{n} u_k * \psi_k)**s \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = x**s
        end function ker
    end function pw_r
    
    !--
    pure function pw_i(this,a,s) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    integer(ip),intent(in)::            s
    real(rp),dimension(0:ubound(a,1)):: o
        !solve \int (\sum_{k=0}^{n} u_k * \psi_k)**s \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = x**s
        end function ker
    end function pw_i
    
    !--
    pure function ex(this,a) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),dimension(0:ubound(a,1)):: o
        !solve \int e^(\sum_{k=0}^{n} c_k*\psi_k) \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = exp(x)
        end function ker
    end function ex
    
    !--
    pure function ln(this,a) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),dimension(0:ubound(a,1)):: o
        !solve \int log(\sum_{k=0}^{n} c_k*\psi_k) \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = log(x)
        end function ker
    end function ln
    
    !--
    pure function dc_func(this,a,lowerfunc,upperfunc,dispt) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    procedure(absf1)::                  lowerfunc,upperfunc
    real(rp),intent(in)::               dispt
    real(rp),dimension(0:ubound(a,1)):: o
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            if(x<dispt) then
                ker = lowerfunc(x)
            else
                ker = upperfunc(x)
            endif
        end function ker
    end function dc_func

    !--
    pure function dc_val(this,a,lower,upper,dispt) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    real(rp),intent(in)::               lower,upper,dispt
    real(rp),dimension(0:ubound(a,1)):: o
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            if(x<dispt) then
                ker = lower
            else
                ker = upper
            endif
        end function ker
    end function dc_val
    
    !<2014-Exploring emerging manycore architectures for 
    !uncertainty quantification through embedded stochastic Galerkin methods>
    !deri = df^i/du_j which is a matrix and lead to [df^i = matmul(df^i/du_j,du_j)]
    !df/du = \sum a_k \phi_k, df^i/du_j = \sum a_k <\phi_i \phi_j \phi_k>
    !for example | df^i/du_j = this%deri(this%op(a_u,df/du)) where u = \sum a_{ui} \phi_i
    pure function deri(this,a_dfdu)
    class(polynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a_dfdu
    real(rp),dimension(0:ubound(a_dfdu,1),0:ubound(a_dfdu,1)):: deri
    integer(ip)::                           i,j
        do j=0,this%truncOd_
            do i=0,this%truncOd_
                deri(i,j) = sum(a_dfdu*this%triMat_(i,j,:))
            enddo
        enddo
    end function deri
    
    !if we only concern the df^i = matmul(df^i/du_j,du_j), it is not necessary to compute df^i/du_j
    !the df/du = \sum a_i \phi_i is enough
    pure function diff(this,a_dfdu,a_du)
    class(polynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a_dfdu,a_du
    real(rp),dimension(0:ubound(a_dfdu,1))::diff
    integer(ip)::                           i,j
        diff = 0._rp
        do i=0,this%truncOd_
            do j=0,this%truncOd_
                diff(i) = diff(i) + sum(a_dfdu*this%triMat_(i,j,:))*a_du(j)
            enddo
        enddo
    end function diff
    
    !----------------------------
    pure integer(ip) function truncOd(this)
    class(polynomialSpace),intent(in):: this
        truncOd = this%truncOd_
    end function truncOd
    
    pure integer(ip) function quadnp(this)
    class(polynomialSpace),intent(in):: this
        quadnp = this%quadnp_
    end function quadnp
    
    pure real(rp) function quadx(this,i)
    class(polynomialSpace),intent(in):: this
    integer(ip),intent(in)::            i
        quadx = this%quadx_(i)
    end function quadx
    
    function triMat(this)
    class(polynomialSpace),target,intent(in)::this
    real(rp),dimension(:,:,:),pointer:: triMat
        triMat => this%triMat_
    end function triMat
    
    !====Test
    pure function op1_ls(this,a,func) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a
    procedure(absf1)::                  func
    real(rp),dimension(0:ubound(a,1)):: o
    real(rp),dimension(:,:),allocatable::mat
    real(rp),dimension(:),allocatable:: quadKernel
    integer(ip)::                       i,j,n,np
    real(rp)::                          x
    
        n = this%truncOd_
        np = 2*this%quadnp_
        allocate(mat(np,0:n), quadKernel(np))
        do j=0,n
            do i=1,np
                x = (i-1)*2._rp/(np-1) - 1._rp
                mat(i,j) = this%quadval(j, x)
            enddo
        enddo
        do i=1,np
            x = (i-1)*2._rp/(np-1) - 1._rp
            quadkernel(i) = func(this%quadval(a, x))
        enddo
        call solvelinearleastsquare(mat, quadkernel)
        o = quadkernel(1:n+1)
        
    end function op1_ls
    
    !--
    pure function op2_ls(this,a,b,func) result(o)
    class(polynomialSpace),intent(in):: this
    real(rp),dimension(0:),intent(in):: a,b
    procedure(absf2)::                  func
    real(rp),dimension(0:ubound(a,1)):: o
    real(rp),dimension(:,:),allocatable::mat
    real(rp),dimension(:),allocatable:: quadKernel
    integer(ip)::                       i,j,n,np
    real(rp)::                          x
    
        n = this%truncOd_
        np = 2*this%quadnp_
        allocate(mat(np,0:n), quadKernel(np))
        do j=0,n
            do i=1,np
                x = (i-1)*2._rp/(np-1) - 1._rp
                mat(i,j) = this%quadval(j, x)
            enddo
        enddo
        do i=1,np
            x = (i-1)*2._rp/(np-1) - 1._rp
            quadkernel(i) = func(this%quadval(a, x), this%quadval(b, x))
        enddo
        call solvelinearleastsquare(mat, quadkernel)
        o = quadkernel(1:n+1)
        
    end function op2_ls    
    
end module polynomialSpace_