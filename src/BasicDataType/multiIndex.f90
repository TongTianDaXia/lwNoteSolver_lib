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
        
    contains
        
        generic::                   init => init_0, init_pvd
        procedure::                 init_0
        procedure::                 init_pvd
        
        procedure::                 vardim
        procedure::                 p
        
        procedure(id),deferred::    idxdim
        
    end type multiIndex
    
    !----------------------------------------------
    !multi_i: all of alpha s.t.{|alpha| = i, i<=p, i \in [0,p]} | idxdim = comb(i+d-1,i)
    type,extends(multiIndex):: multiTriIndex
    
        private
        integer(ip)::               h = 0, t = 0
        logical(lp)::               more = .false.
        
    contains
        procedure::                 idxdim => idxdim_tri
        procedure::                 traverse => traverse_tri
    end type multiTriIndex

    !multi_i: all of alpha s.t.{|alpha|_0 = i, i<=|alpha|<=p, i range [0,d]} 
    !idxdim = comb(p+d,p) - (comb(i+d,i)-comb(i+d-1,i)) for i>0
    !idxdim = 1 for i=0
    !refer to <Global sensitivity analysis using polynomial chaos expansions>
    type,extends(multiIndex):: multiAnovaIndex
    
        private
        integer(ip)::               cp = -1
        integer(ip)::               h = 0, t = 0
        logical(lp)::               combmore = .false.
        logical(lp)::               compmore = .false.
        integer(ip),dimension(:),&
        allocatable::               iv,a
        
    contains
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
    integer(ip)::                           p,d
        
        p = this%p_; d = this%vardim_
        if(i/=this%i_) then
            call this%init(p,d) !can't put this%p_, init first let p=minip
            this%i_ = i
        end if
        call compositionNext(p,d,alpha,this%more,this%h,this%t)
        last = .not.this%more
        if(last) call this%init(p,d)
        
    end subroutine traverse_tri
    
    !-----------------------------------------------------------------    
    pure integer(ip) function idxdim_Anova(this,i) result(d)
    class(multiAnovaIndex),intent(in)::     this
    integer(ip),intent(in)::                i
        d = merge(1_ip, nint(binomialCoef(this%p_+this%vardim_,this%p_)) - nint(binomialCoef(i+d-1,i-1)), i==0)
    end function idxdim_Anova
    
    pure subroutine traverse_Anova(this,i,alpha,last)
    class(multiAnovaIndex),intent(inout)::      this
    integer(ip),intent(in)::                    i
    integer(ip),dimension(:),intent(out)::      alpha
    logical(lp),intent(out)::                   last
    integer(ip)::                               p,d,j
        
        p = this%p_; d = this%vardim_
        last = .false.
        alpha = 0
        
        if(i==0) then
        
            alpha = 0
            last = .true.
            
        else
        
            if(i/=this%i_) then
                call this%init(p,d)
                this%i_ = i
                if(allocated(this%iv)) deallocate(this%iv)
                if(allocated(this%a)) deallocate(this%a)
                allocate(this%iv(i),this%a(i))
            end if
        
            if(p<i) then
                alpha = -1
                last = .true.
            elseif(p==i) then
                call combinationNext(this%vardim_,i,this%iv,this%combmore)
                do j=1,i
                    alpha(this%iv(j)) = 1
                enddo
                last = .not. this%combmore
            elseif(p>i) then
                !traverse all of possible (vardim,i)
                !if no more composition and no more cp, then iterate combination
                if(.not.this%compmore .and. this%cp==-1) &
                    call combinationNext(this%vardim_, i, this%iv, this%combmore)
                
                !under this idxvector, traverse all of possible element
                if(this%cp==-1) this%cp = 0 !init
                call compositionNext(this%cp, i, this%a, this%compmore, this%h, this%t)
                !composition make sulre |this%a|<=p-i, then scaled to i<=|this%a|<=p
                do j=1,i
                    alpha(this%iv(j)) = this%a(j) + 1
                enddo
                
                if(.not.this%compmore) this%cp = this%cp + 1    !if last under this cp, then plus cp by 1
                if(this%cp > p-i) this%cp = -1  !if cp over p-i, then ruin it
                last = .not.this%combmore .and. .not.this%compmore .and. this%cp==-1
            endif
        endif
        if(last) call this%init(p,d)
        
    end subroutine traverse_Anova
    
end module multiIndex_