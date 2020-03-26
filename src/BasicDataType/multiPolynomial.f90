module multiPolynomial_
use constants
use polynomial_
implicit none

	private
	public:: multiPolynomial

	!--
	type:: multiPolynomial
	
		type(polynomial),dimension(:),allocatable::	po
		
	contains
	
		generic::			init => init_dim
		procedure::			init_dim
	
	end type multiPolynomial
	
	
contains


	!--
	subroutine init_dim(this, dim)
	class(multiPolynomial),intent(out)::	this
	integer(ip),intent(in)::				dim
	
		allocate(this%po(dim))
	
	end subroutine init_dim

end module multiPolynomial_