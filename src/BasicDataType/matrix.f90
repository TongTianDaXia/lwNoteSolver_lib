module matrix_
use constants
use vector_
implicit none

    private
    
    public:: matrix      
    public:: operator(*), operator(+), operator(-)
    
    !--
    type:: matrix
	
        private
        real(rp),allocatable,dimension(:,:)::   mat
		
    contains
	
		!--
        generic::           init => matrix_zero, matrix_2rank
        procedure::         matrix_zero
        procedure::         matrix_2rank
		
		!--
        procedure::         init_movealloc
        
        !--
        procedure::         allocstat
        
        !sval
        generic::           sval => sval_columnVector, sval_locVal, sval_RowAR
        procedure::         sval_columnVector
        procedure::         sval_RowAR
        procedure::         sval_locVal
        
        !member function
        procedure::         ptr	=> getMatptr
        procedure::         getMatptr
		
        procedure::         vals => getMatvals
        procedure::         getMatvals
        
        procedure::         Cptr    => getColumnptr
        procedure::         getColumnptr
		
        procedure::         Rptr    => getRowptr
        procedure::         getRowptr
		
        !----
        !for now i don't want offering the function to get vals of Column and Row
        procedure::         vector => getVector
        procedure::         getVector
		
        procedure::         val	=> getval
        procedure::         getVal
		
        procedure::         cdim => getcolumndim
        procedure::         getColumnDim
		
        procedure::         rdim    => getRowdim
        procedure::         getRowDim
        
        !-special member func
        procedure::         diag
        procedure::         trace
		
    end type matrix

    !--------
    interface operator(+)
        procedure::  vmplus
        procedure::  mvplus
    end interface
    
    interface operator(-)
        procedure::  negm
        procedure::  vmminus
        procedure::  mvminus
    end interface
    
    interface operator(*)
        procedure::  smproduct
        procedure::  msproduct
        procedure::  vmproduct
        procedure::  mvproduct
        procedure::  mmproduct
    end interface
    
!--
contains

	!--
    elemental subroutine matrix_zero(this,cd,rd)
    class(matrix),intent(out)::			this
    integer(ip),intent(in)::			cd
    integer(ip),intent(in),optional::	rd
    integer(ip)::						cd_,rd_
	
        cd_ = cd
        rd_ = merge(rd, cd, present(rd))
        allocate(this%mat(cd_, rd_))
        this%mat = zero
		
    end subroutine matrix_zero
    
    !--
    pure subroutine matrix_2rank(this,vals)
    class(matrix),intent(out)::			this
    real(rp),dimension(:,:),intent(in)::vals
	
        allocate(this%mat, source = vals)
		
    end subroutine matrix_2rank
    
    !--
    pure subroutine init_MoveAlloc(this,vals)
    class(matrix),intent(out)::         this
    real(rp),allocatable,dimension(:,:),&
    intent(inout)::						vals
	
        call move_alloc(vals, this%mat)
		
    end subroutine init_MoveAlloc
    
    
    logical(lp) function allocstat(this)
    class(matrix),intent(in)::          this
	
        allocstat = allocated(this%mat)
		
    end function allocstat
    
    
	!---member function 
    impure elemental subroutine sval_columnVector(this,Column_i,vct)
    class(matrix),intent(inout)::       this
    class(vector),intent(in)::          vct
    integer(ip),intent(in)::            Column_i
	
        this%mat(:, Column_i) = vct%ptr()
		
    end subroutine sval_columnVector
	
    !--
    pure subroutine sval_RowAR(this,Row_j,ar)
    class(matrix),intent(inout)::       this
    integer(ip),intent(in)::            Row_j
    real(rp),dimension(:),intent(in)::  ar
	
        this%mat(Row_j, :) = ar
		
    end subroutine sval_RowAR
    
	!--
    pure subroutine sval_locVal(this,Column_i,Row_j,val)
    class(matrix),intent(inout)::       this
    integer(ip),intent(in)::            Column_i,Row_j
    real(rp),intent(in)::               val
	
        this%mat(Column_i, Row_j) = val
		
    end subroutine sval_locVal
   
    !--
    function getMatptr(this) result(mat)
    class(matrix),target,intent(in)::   this
    real(rp),pointer,dimension(:,:)::   mat
	
        mat => this%mat
		
    end function getMatptr
    
    !--
    pure function getMatvals(this) result(mat)
    class(matrix),intent(in)::          this
    real(rp),allocatable,dimension(:,:)::mat
	
        allocate(mat, source = this%mat)
		
    end function getMatvals
   
    
    !--
    function getColumnptr(this,Row_j) result(ar)
    class(matrix),target,intent(in)::   this
    integer(ip),intent(in)::            Row_j
    real(rp),pointer,dimension(:)::     ar
	
        ar => this%mat(:, Row_j)
		
    end function getColumnptr
   
    !--
    function getRowptr(this,Column_i) result(ar)
    class(matrix),target,intent(in)::   this
    integer(ip),intent(in)::            Column_i
    real(rp),pointer,dimension(:)::     ar
	
        ar => this%mat(Column_i, :)
		
    end function getRowptr
    
    !--
    elemental type(vector) function getVector(this, Row_j ) result(vt)
    class(matrix),intent(in)::          this
    integer(ip),intent(in)::            Row_j
	
        call vt%init(this%mat(:, Row_j))
		
    end function getVector
    
    !--
    pure real(rp) function getVal(this ,Column_i, Row_j ) result(val)
    class(matrix),intent(in)::          this
    integer(ip),intent(in)::            Column_i,Row_j
	
        val = this%mat(Column_i,Row_j)
		
    end function getVal
    
    !--
    pure integer(ip) function getColumnDim(this) result(d)
    class(matrix),intent(in)::          this
	
        d = size(this%mat, 1)
		
    end function getColumnDim
    
    !--
    pure integer(ip) function getRowDim(this) result(d)
    class(matrix),intent(in)::  this
	
        d = size(this%mat, 2)
		
    end function getRowDim
    
    
	!--
    elemental type(vector) function diag(this)
    class(matrix),intent(in)::  this
    integer(ip)::               d,i
	
        d = min(this%cdim(), this%rdim())
        call diag%init(d)
        do i=1,d
            call diag%sval(i, this%mat(i,i))
        end do
		
    end function diag
    
	!--
    elemental real(rp) function trace(this)
    class(matrix),intent(in)::  this
    integer(ip)::               d,i
	
        d = min(this%cdim(),this%rdim())
        trace = zero
        do i=1,d
            trace = trace + this%mat(i,i)
        enddo
		
    end function trace
    
	!--------------------------------------------
    !--op(+)
    elemental type(matrix) function mvplus(lhs,rhs) result(p)
    type(matrix),intent(in)::   lhs
    type(vector),intent(in)::   rhs
    integer(ip)::               d,i,j
        d = min(rhs%dim(), lhs%cdim(), lhs%rdim())
        allocate(p%mat, source=lhs%mat)
        forall(i=1:d,j=1:d,i==j) p%mat(i,j) = p%mat(i,j) + rhs%val(i)
    end function mvplus
    
    elemental type(matrix) function vmplus(lhs,rhs) result(p)
    type(matrix),intent(in)::   rhs
    type(vector),intent(in)::   lhs
        p = rhs + lhs
    end function vmplus
    
    
    !--op(-)
    elemental type(matrix) function negm(rhs) result(m)
    type(matrix),intent(in)::   rhs
        m%mat = -rhs%mat
    end function negm
    
	!--
    elemental type(matrix) function mvminus(lhs,rhs) result(p)
    type(matrix),intent(in)::   lhs
    type(vector),intent(in)::   rhs
        p = lhs + (-rhs)
    end function mvminus
    
	!--
    elemental type(matrix) function vmminus(lhs,rhs) result(p)
    type(matrix),intent(in)::   rhs
    type(vector),intent(in)::   lhs
        p = - (rhs - lhs)
    end function vmminus
    
    !--
    elemental type(matrix) function smproduct(lhs,rhs) result(p)
    real(rp),intent(in)::       lhs
    type(matrix),intent(in)::   rhs
        allocate(p%mat,source = rhs%mat * lhs)
    end function smproduct
    
    !--
    elemental type(matrix) function msproduct(lhs,rhs) result(p)
    type(matrix),intent(in)::   lhs
    real(rp),intent(in)::       rhs
        allocate(p%mat,source = lhs%mat * rhs)
    end function msproduct
    
    !--
    impure elemental type(vector) function mvproduct(lhs,rhs) result(p)
    type(matrix),intent(in)::   lhs
    type(vector),intent(in)::   rhs
        call p%init(matmul(lhs%mat, rhs%ptr()))
    end function mvproduct
    
    !--
    impure elemental type(vector) function vmproduct(lhs,rhs) result(p)
    type(vector),intent(in)::   lhs
    type(matrix),intent(in)::   rhs
    real(rp),dimension(:,:),allocatable:: v
        allocate(v(1,lhs%dim()))
        v(1,:) = lhs%ptr()
        v = matmul(v,rhs%mat)
        call p%init(v(1,:))
    end function vmproduct
    
    !--
    elemental type(matrix) function mmproduct(lhs,rhs) result(p)
    type(matrix),intent(in)::   lhs,rhs
        allocate(p%mat,source = matmul(lhs%mat,rhs%mat))
    end function mmproduct

end module matrix_