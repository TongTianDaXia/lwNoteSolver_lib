!tensor construction: |alpha|_1<=p; total construction: |alpha|_\infty <= p; hyperbolic construction: \prod{\alpha_i+1} <= p+1
!more info refer to <不确定量化的高精度数值方法和理论>
!for now, only tensor construction is under consideration
module multiIndex_
use constants
use arrayOpsLib
use SpecialFunctionLib
implicit none
    
    private
    public::    multiIndex
    public::    multiTriIndex
    public::    multiAnovaIndex
    
    !alpha is the element of multiIndex | e.g. alpha=(1,0,0,2)
    !dim(multiIndex_i) = idxdim; the number of alpha belong to multiIndex_i
    !dim(alpha) = vardim, the number of variate
    !alpha_i range [0,p]
    type,abstract:: multiIndex
        
        private
        integer(ip)::               p_ = minip
        integer(ip)::               vardim_ = minip

        !for traverse
        integer(ip)::               i_ = minip
        integer(ip)::               h = 0, t = 0
        logical(lp)::               more = .false.
        
    contains
        
        generic::                   init => init_0, init_pvd
        procedure::                 init_0
        procedure::                 init_pvd
        
        procedure::                 vardim
        procedure::                 p
        
        procedure(id),deferred::idxdim
        
    end type multiIndex
    
    !----------------------------------------------
    !multi_i: all of alpha s.t.{|alpha| = i, i<=p, i \in [0,p]} | idxdim = comb(i+d-1,i)
    type,extends(multiIndex):: multiTriIndex
        private
    contains
        procedure::                 idxdim => idxdim_tri
        procedure::                 traverse => traverse_tri
    end type multiTriIndex

    !multi_i: all of alpha s.t.{|alpha|_0 = i, i<=|alpha|<=p, i range [0,d]} 
    !idxdim = comb(p+d,p) - (comb(i+d,i)-comb(i+d-1,i)) for i>0
    !idxdim = 1 for i=0
    type,extends(multiIndex):: multiAnovaIndex
        private
        integer(ip),dimension(:),allocatable:: iv
    contains
        procedure::                 init_pvd => init_pvd_anova
        procedure::                 idxdim => idxdim_Anova
        procedure::                 traverse => traverse_Anova
    end type multiAnovaIndex
    
    
    !---------
    abstract interface
        pure integer(ip) function id(this,i) result(d)
        import:: ip,multiIndex
        class(multiIndex),intent(in)::this
        integer(ip),intent(in)::      i
        end function id
    end interface
    
    
contains
    
    
    pure subroutine init_0(this)
    class(multiIndex),intent(out)::     this
    end subroutine init_0
    
    pure subroutine init_pvd(this,p,vardim)
    class(multiIndex),intent(out)::     this
    integer(ip),intent(in)::            p,vardim
        this%p_ = p
        this%vardim_ = vardim
    end subroutine init_pvd
    
    pure integer(ip) function vardim(this)
    class(multiIndex),intent(in)::      this
        vardim = this%vardim_
    end function vardim
    
    pure integer(ip) function p(this)
    class(multiIndex),intent(in)::      this
        p = this%p_
    end function p
    
    !---------------------------------------------------------------
    pure integer(ip) function idxdim_tri(this,i) result(d)
    class(multiTriIndex),intent(in)::       this
    integer(ip),intent(in)::                i
        d = nint(BinomialCoef(i+this%vardim_-1, i))
    end function idxdim_tri
    
    pure subroutine traverse_tri(this,i,alpha,last)
    class(multiTriIndex),intent(inout)::    this
    integer(ip),intent(in)::                i
    integer(ip),dimension(:),intent(inout)::alpha
    logical(lp),intent(out)::               last
        
        if(i/=this%i_) then
            call this%init(this%p_,this%vardim_)
            this%i_ = i
        end if
        call compositionNext(this%p_,this%vardim_,alpha,this%more,this%h,this%t)
        last = .not.this%more
        if(last) call this%init(this%p_,this%vardim_)
        
    end subroutine traverse_tri
    
    !-----------------------------------------------------------------
    pure subroutine init_pvd_anova(this,p,vardim)
    class(multiAnovaIndex),intent(out)::    this
    integer(ip),intent(in)::                p,vardim
        this%p_ = p
        this%vardim_ = vardim
        allocate(this%iv(vardim)); this%iv = minip
    end subroutine init_pvd_anova
    
    pure integer(ip) function idxdim_Anova(this,i) result(d)
    class(multiAnovaIndex),intent(in)::     this
    integer(ip),intent(in)::                i
        d = merge(1_ip, nint(binomialCoef(this%p_+this%vardim_,this%p_)) - nint(binomialCoef(i+d-1,i-1)), i==0)
    end function idxdim_Anova
    
    pure subroutine traverse_Anova(this,i,alpha,last)
    class(multiAnovaIndex),intent(inout)::      this
    integer(ip),intent(in)::                    i
    integer(ip),dimension(:),intent(inout)::    alpha
    logical(lp),intent(out)::                   last
        
        last = .false.
        if(i==0) then
            call this%init(this%p_,this%vardim_)
            alpha = 0
            last = .true.
        else
            if(i/=this%i_) then
                call this%init(this%p_,this%vardim_)
                this%i_ = i
            end if
            call compositionNext(i,this%vardim_,this%iv
            ,this%more,this%h,this%t)    
        endif
        
        
        !call compositionNext(this%p_,this%vardim_,alpha,this%more,this%h,this%t)
        !last = .not.this%more
        !if(last) call this%init(this%p_,this%vardim_)
        
    end subroutine traverse_Anova
    
end module multiIndex_