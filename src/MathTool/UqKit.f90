module UqKit_
use constants
use samplinglib
use IntegrationLib
use multiIndex_
use polynomial_
use specialPolynomial_
implicit none

	private
	public:: UqKit


	type:: UqKit

		private

		character(cl)::				tp_

		!quadature
		real(rp),dimension(:,:),&
		allocatable::				x_
		real(rp),dimension(:),&
		allocatable::				w_

		!polynomial expansion
		integer(ip),dimension(:,:),&
		allocatable::				alpha_
		real(rp),dimension(:),&
		allocatable::				coefpo_
		type(polynomial),dimension(:),&
		allocatable::				bas_


	contains

		generic::		init => init_tp
		procedure::		init_tp

		!--
		procedure::		makeSparseGrid
		procedure::		integrate => integrate_uk
		! procedure::		expct
		! procedure::		std

		!--
		procedure::		makeCoefpo
		! procedure::		makepdf

		!--
		procedure::		dim
		procedure::		nqp	!number of quadrature points
		procedure::		npo	!number of polynomial expansion

	end type UqKit

	!--
	abstract interface
		real(rp) function fc(x)
		import:: rp
		real(rp),dimension(:),intent(in):: 	x
		end function fc
	end interface

contains

	!--
	subroutine init_tp(this, tp)
	class(UqKit),intent(out)::	this
	character(*),intent(in)::	tp

		this%tp_ = tp

	end subroutine init_tp

	!--
	subroutine makeSparseGrid(this, dim, lvl, quadrule)
	class(UqKit),intent(inout)::	this
	integer(ip),intent(in)::		dim, lvl
	character(*),intent(in)::		quadrule
	integer(ip)::					n

		this%tp_ = quadrule

		n = sparseGridSize(dim, lvl, quadrule)
		allocate(this%x_(dim, n), this%w_(n))
		call sparseGrid(lvl, quadrule, this%x_, this%w_)

		!--
		if(quadrule=='gh') then

			!in this condition, native (x, w) is the quadrule for \int f(x) e^{-x^2} dx
			!convert to e^{-x^2}, (x, w) -> (sqrt(2)x, sqrt(2)^dim w), which just the quadrule for 'prob'
			!convert to \frac{1}{\sqrt(2 \pi)^dim} e^{-x^2}, (x, w) -> (sqrt(2)x, w/\sqrt(\pi)^dim)
			this%x_ = sqrt(2._rp)*this%x_
			this%w_ = this%w_/srpi**dim

		elseif(quadrule=='gl' .or. quadrule=='cc') then

			!in this condition, native (x, w) is the quadrule for \int_{-1}^{1} f(x) dx
			!convert to probability \frac{1}{2}, (x, w) -> (x, w/2^dim)
			this%w_ = this%w_/2._rp**dim

		endif

	end subroutine makeSparseGrid

	!--
	real(rp) function integrate_uk(this, f, mu, sigma) result(itg)
	class(UqKit),intent(in)::			this
	procedure(fc)::						f
	real(rp),dimension(:),optional,&
	intent(in)::						mu, sigma
	integer(ip)::						i
	real(rp),dimension(:),allocatable::	mu0, sigma0

		allocate(mu0(this%dim()), sigma0(this%dim()))
		mu0 = merge(mu, 0._rp, present(mu))
		sigma0 = merge(sigma, 1._rp, present(sigma))

		itg = 0._rp
		do i=1,this%nqp()
			itg = itg + f(sigma0*this%x_(:, i) + mu0)*this%w_(i)
		end do

	end function integrate_uk

	!--
	! real(rp) function expec()

	!--
	subroutine makeCoefpo(this, f, maxp, maxcovdim)
	class(UqKit),intent(inout)::		this
	procedure(fc)::						f
	integer(ip),intent(in)::			maxp, maxcovdim
	integer(ip)::						i, j, npo, dim
	type(multiAnovaIndex)::				AnovaIdx
	logical(lp)::						more

		dim = this%dim()
		if(maxcovdim > dim .or. maxcovdim < 0) &
		stop 'error: UqKit%makeCoefpo get an error maxcovdim'

		!make alpha
		call AnovaIdx%init(maxp, dim)
		npo = 0
		do i=0,maxcovdim
			npo = npo + AnovaIdx%idxdim(i)
		end do

		allocate(this%alpha_(this%dim(), npo), this%coefpo_(npo), this%bas_(maxp))
		npo = 0
		do i=0,maxcovdim
			more = .true.
			do while(more)
				npo = npo + 1
				call AnovaIdx%traverse(i, this%alpha_(:, npo), more)
			enddo
		enddo

		!make stardard basis polynomial
		if(this%tp_ == 'gh') then
			this%bas_ = normalHermitePolynomialSet(maxp, 'prob')
			this%bas_ = sqrt(sqrt(2._rp)*srpi)*this%bas_
		else
			stop 'error: Uqkit~'
		endif

		!make coef
		this%coefpo_ = 0._rp
		do i=1,npo
			this%coefpo_(i) = this%integrate(projker)
		enddo

	contains

		real(rp) function projker(x)
		real(rp),dimension(:),intent(in)::	x

			projker = f(x)*multiPolynominal(this%bas_, this%alpha_(:,i), x)

		end function projker

	end subroutine makeCoefpo


	!--
	pure integer(ip) function dim(this)
	class(UqKit),intent(in):: 	this

		dim = size(this%x_, 1)

	end function dim

	!--
	pure integer(ip) function nqp(this)
	class(UqKit),intent(in):: 	this

		nqp = size(this%x_, 2)

	end function nqp

	!--
	pure integer(ip) function npo(this)
	class(UqKit),intent(in):: 	this

		npo = size(this%coefpo_)

	end function npo

end module UqKit_
