module vectorvector_
use constants
use vector_
implicit none

	private
	public:: vectorvector
	
	type:: vectorvector
	
		private
		type(vector),dimension(:),allocatable::	vc
	
	contains
	
		!--
        generic::       init    =>  init_n,     &
                                    init_size
        procedure::     init_n
        procedure::     init_size
		
		!--
        procedure::     init_movealloc
        
        !--
        generic::       sval    => sval_val, sval_vals
        procedure::     sval_val
        procedure::     sval_vals
        
        !---
        procedure::     val => val_ij
        procedure::     val_ij
		
		!--
        generic::       ptr => ptr_i
        procedure::     ptr_i
		
		!--
		generic::		vsize => jsize, size_i
		procedure::		jsize
		procedure::		size_i
	
	end type vectorvector
	
!--
contains

	!--
	subroutine init_n(this,vsize,nvc)
    class(vectorvector),intent(out)::		this
    integer(ip),intent(in)::				vsize,nvc
    integer(ip)::							i
	
        allocate(this%vc(nvc))
        do i = 1,nvc
            call this%vc(i)%init(vsize)
        enddo
		
    end subroutine init_n
	
	!--
    subroutine init_size(this,sizes)
    class(vectorvector),intent(out)::		this
    integer(ip),dimension(:),intent(in)::	sizes
    integer(ip)::							i
	
        allocate(this%vc(size(sizes)))
        do i = 1,size(sizes)
            call this%vc(i)%init(sizes(i))
        enddo
		
    end subroutine init_size
	
	!--
    subroutine init_movealloc(this,vcvc)
    class(vectorvector),intent(out)::   this
    class(vectorvector),intent(inout):: vcvc
    integer(ip)::						i,n
	
        n = size(vcvc%vc)
        allocate(this%vc(n))
        do i=1,n
            call this%vc(i)%init_movealloc(vcvc%vc(i))
        enddo
		
    end subroutine init_movealloc

    !--
    elemental subroutine sval_val(this, i, j, val)
    class(vectorvector),intent(inout):: this
    integer(ip),intent(in)::			i, j
	real(rp),intent(in)::				val
	
        call this%vc(j)%sval(i, val)
		
    end subroutine sval_val
	
	!--
    pure subroutine sval_vals(this, j, vals)
    class(vectorvector),intent(inout)::	this
    integer(ip),intent(in)::			j
    real(rp),dimension(:),intent(in)::	vals
	
		call this%vc(j)%sval(vals)
		
    end subroutine sval_vals
    
    !--
    elemental function val_ij(this, i, j)
    class(vectorvector),intent(in)::	this
    integer(ip),intent(in)::			i,j
    integer(ip)::						val_ij
	
        val_ij = this%vc(j)%val(i)
		
    end function val_ij
    
    !--
    function ptr_i(this,i) result(ptr)
    class(vectorvector),target,intent(in)::	this
    integer(ip),intent(in)::				i
    type(vector),pointer::					ptr
	
        ptr => this%vc(i)
		
    end function ptr_i
	
	!--
    elemental function jsize(this) result(s)
    class(vectorvector),intent(in)::	this
    integer(ip)::						s
	
        s = size(this%vc)
		
    end function jsize
    
    !--
    elemental function size_i(this,i) result(s)
    class(vectorvector),intent(in)::    this
    integer(ip),intent(in)::            i
    integer(ip)::                       s
	
        s =  this%vc(i)%dim()
		
    end function size_i
	
end module vectorvector_