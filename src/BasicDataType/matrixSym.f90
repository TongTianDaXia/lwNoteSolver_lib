module matrixSym_
use constants
use arrayOpsLib
implicit none

	private
	public:: matrixSym
	
	!sz*sz symmetric matrix
	type:: matrixSym

		private
		integer(ip)::			sz_
		
		real(rp),dimension(:),&
		allocatable::			mat_
		
	contains
		
		procedure::		init
		procedure::		sval
		
		!--
		procedure::		val => val_ij
		procedure::		vals => val_mat
		procedure::		mat => val_mat
		
		!--
		procedure::		idx
	
	end type
	
	
contains

	elemental subroutine init(this, sz)
	class(matrixSym),intent(out)::		this
	integer(ip),intent(in)::			sz
	
		this%sz_ = sz
		allocate(this%mat_(pidb2(sz*(sz+1))))

	end subroutine init
	
	!--
	elemental subroutine sval(this, i, j, val)
	class(matrixSym),intent(inout)::	this
	integer(ip),intent(in)::			i, j
	real(rp),intent(in)::				val
	
		this%mat_(this%idx(i, j)) = val
	
	end subroutine sval
	
	elemental real(rp) function val_ij(this, i, j)
	class(matrixSym),intent(in)::	this
	integer(ip),intent(in)::		i, j
	
		val_ij = this%mat_(this%idx(i, j))
	
	end function val_ij
	
	!--
	function val_mat(this) result(val)
	class(matrixSym),intent(in)::	this
	integer(ip)::					i, j
	real(rp),dimension(:,:),&
	allocatable::					val
	
		allocate(val(this%sz_, this%sz_))
	
		do j=1,this%sz_
			do i=1,this%sz_
				val(i, j) = this%mat_(this%idx(i, j))
			enddo
		enddo
	
	end function val_mat
	
	!top right
	elemental integer(ip) function idx(this, i, j)
	class(matrixSym),intent(in)::	this
	integer(ip),intent(in)::		i, j
	integer(ip)::					ii, jj
	
		if(i <= j) then
			ii = i
			jj = j
		else
			ii = j
			jj = i
		endif
		
		idx = pidb2(jj*(jj - 1)) + ii
		
	end function idx
	
end module matrixSym_