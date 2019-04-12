module arrayOpsLib
use constants, only:ip,rp,lp,zero,swap
implicit none

    private
    
    !--operator for vector and matrix
    public:: operator(.ip.),operator(.op.),operator(.cpv.),operator(.cps.)
    public:: operator(-),operator(+),operator(*),operator(.eqv.),operator(.eqvs.)
    
    !--operation based on the array
    public:: magSqr, mag, angle, normalize, para, orth, rot2, rot2deg, norm, polyval
    public:: trace, diag, unitary, sort, cumprod, trans, repmat, bondcol, bondrow
    
    !--indice ops
    public:: compositionNext,combinationNext,colexNext
    
    !--compress stored row, for sparse matrix
    public:: diagIxCsr
    
    
    !==============================================
    interface operator(*)
        procedure:: mvinnerproduct
        procedure:: vminnerproduct
    end interface 
    
    interface operator(.ip.)
        procedure:: vvinnerproduct
        procedure:: mminnerproduct
        procedure:: mvinnerproduct
        procedure:: vminnerproduct
    end interface
    
    interface operator(.op.)
        procedure:: vvouterproduct
    end interface
    
    interface operator(.cpv.)
        procedure:: vvcrossproduct3
    end interface
    
    interface operator(.cps.)
        procedure:: vvcrossproduct2
    end interface
    
    interface operator(+)
        procedure:: mvplus
        procedure:: vmplus
    end interface

    interface operator(-)
        procedure:: mvminus
        procedure:: vmminus
    end interface
    
    interface operator(.eqv.)
        procedure:: eqv_r
        procedure:: eqv_i
    end interface
    
    interface operator(.eqvs.)
        procedure:: eqvs
    end interface
    
    
    !---------------
    interface polyval
        procedure:: polyval_horner
        procedure:: polyval_mat_horner
        procedure:: polyval_matv
    end interface polyval
    
    interface rot2
        procedure:: rot2Cs
        procedure:: rot2Theta
    end interface rot2
    
    interface diag
        procedure:: diagCreateMatrix
        procedure:: diagExtractElement
    end interface diag
    
    interface sort  !bubble sort
        procedure:: sort1_i
        procedure:: sort2_i
        procedure:: sort1_r
        procedure:: sort2_r
    end interface sort
    
    interface trans
        procedure:: trans_rvm
        procedure:: trans_ivm
        procedure:: trans_rmm
    end interface trans
    
    interface repmat
        procedure:: repmat_rv
        procedure:: repmat_iv
        procedure:: repmat_rm
        procedure:: repmat_im
    end interface repmat
    
    interface bondcol
        procedure:: bondcol_vvr
        procedure:: bondcol_vvvr
    end interface bondcol
    
    interface bondrow
        procedure:: bondrow_vvr
        procedure:: bondrow_vvvr
    end interface bondrow
    
    interface colexNext
        procedure:: colexNext_isotropy
        procedure:: colexNext_anisotropy
    end interface colexNext
    
!=======================================================
contains
!=======================================================

    pure real(rp) function mag(v)
    real(rp),dimension(:),intent(in)::  v
        mag = norm2(v)
    end function mag
    
    pure real(rp) function magSqr(v) result(ms)
    real(rp),dimension(:),intent(in)::  v
        ms = sum(v**2)
    end function magSqr
    
    pure real(rp) function angle(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
        angle = acos((v1.ip.v2)/mag(v1)/mag(v2))
    end function angle
    
    pure function normalize(v)
    real(rp),dimension(:),intent(in)::  v
    real(rp),dimension(size(v))::       normalize
        normalize = v / mag(v)
    end function normalize
    
    pure function para(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
    real(rp),dimension(size(v1))::      para
        para = (v1.ip.v2)*normalize(v2)
    end function para
    
    pure function orth(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
    real(rp),dimension(size(v1))::      orth
        orth = v1 - para(v1,v2)
    end function orth
    
    !cs means cos and sin, v is the vector
    pure function rot2Cs(v,cs) result(r)
    real(rp),dimension(2),intent(in)::  cs,v
    real(rp),dimension(2)::             r
        r(1) = cs(1)*v(1) - cs(2)*v(2)
        r(2) = cs(2)*v(1) + cs(1)*v(2)
    end function rot2Cs
    
    pure function rot2Theta(v,theta) result(r)
    real(rp),dimension(2),intent(in)::  v
    real(rp),intent(in)::               theta
    real(rp),dimension(2)::             r
        r(1) = cos(theta)*v(1) - sin(theta)*v(2)
        r(2) = sin(theta)*v(1) + cos(theta)*v(2)
    end function rot2Theta
    
    pure function rot2deg(v,theta) result(r)
    real(rp),dimension(2),intent(in)::  v
    real(rp),intent(in)::               theta
    real(rp),dimension(2)::             r
        r(1) = cosd(theta)*v(1) - sind(theta)*v(2)
        r(2) = sind(theta)*v(1) + cosd(theta)*v(2)
    end function rot2deg
    
    !Lp norm
    pure real(rp) function norm(v,p)
    real(rp),dimension(:),intent(in)::  v
    integer(ip),intent(in)::            p
    integer(ip)::                       i
        if(p==2) then
            norm = norm2(v)
        elseif(p==1) then
            norm = sum(abs(v))
        elseif(p==0) then
            norm = count(v /= [(0._rp,i=1,size(v))])
        elseif(p>2) then
            norm = sum(abs(v)**p)**(1._rp/p)
        else
            norm = maxval(abs(v))
        endif
    end function norm
    
    !---------
    !pure real(rp) function polyval(a,x)
    !real(rp),dimension(0:),intent(in):: a
    !real(rp),intent(in)::               x
    !integer(ip)::                       i
    !    polyval = 0._rp
    !    if(abs(x)<1._rp) then
    !        do i=ubound(a,1),0,-1
    !            polyval = polyval + a(i)*x**i
    !        enddo
    !    else
    !        do i=0,ubound(a,1)
    !            polyval = polyval + a(i)*x**i
    !        enddo
    !    endif
    !end function polyval
    !--refer to https://en.wikipedia.org/wiki/Horner%27s_method#cite_note-9
    pure real(rp) function polyval_horner(a,x) result(val)
    real(rp),dimension(0:),intent(in):: a
    real(rp),intent(in)::               x
    integer(ip)::                       i
        val = 0._rp
        do i=ubound(a,1),0,-1
            val = val*x + a(i)
        enddo
    end function polyval_horner
    
    !compute p(x)
    pure function polyval_mat_horner(a,x) result(val)
    real(rp),dimension(0:),intent(in)::         a
    real(rp),dimension(:,:),intent(in)::        x
    real(rp),dimension(size(x,1),size(x,2))::   val
    integer(ip)::                               i,j,n
        n = size(x,1)
        val = 0._rp
        do i=ubound(a,1),0,-1
            val = (val.ip.x) + [(a(i),j=1,n)]
        enddo
    end function polyval_mat_horner
    
    !compute p(x)v
    pure function polyval_matv(a,x,v) result(val)
    real(rp),dimension(0:),intent(in)::     a
    real(rp),dimension(:,:),intent(in)::    x
    real(rp),dimension(:),intent(in)::      v
    real(rp),dimension(size(v))::           val,xv
    integer(ip)::                           i
        xv = v
        val = a(0)*v
        do i=1,ubound(a,1)
            xv = x.ip.xv
            val = val + a(i)*xv
        enddo
    end function polyval_matv
    
    
    !---------
    pure function diagCreateMatrix(m,k)
    real(rp),dimension(:),intent(in)::          m
    integer(ip),intent(in)::                    k
    real(rp),dimension(:,:),allocatable::       diagCreateMatrix
    integer(ip)::                               i,j,n,p
        p = abs(k); n = size(m)
        allocate(diagCreateMatrix(p+n,p+n))
        diagCreateMatrix = 0._rp
        if(k>=0) then
            do i=1,n
                diagCreateMatrix(i,i+k) = m(i)
            enddo
        elseif(k<0) then
            do j=1,n
                diagCreateMatrix(j+p,j) = m(j)
            enddo
        endif 
    end function diagCreateMatrix
    
    !--
    pure function diagExtractElement(m)
    real(rp),dimension(:,:),intent(in)::    m
    real(rp),dimension(min(size(m,dim=1),&
    size(m,dim=2)))::                       diagExtractElement
    integer(ip)::                           i
        forall(i=1:size(diagExtractElement)) diagExtractElement(i)=m(i,i)
    end function diagExtractElement
    
    
    !--create unitary matrix
    pure function unitary(n)
    integer(ip),intent(in)::    n
    real(rp),dimension(n,n)::   unitary
    integer(ip)::               i
        unitary = zero
        do i=1,n
            unitary(i,i) = 1._rp
        enddo
    end function unitary
    
    !----
    pure real(rp) function trace(m)
    real(rp),dimension(:,:),intent(in)::    m
    integer(ip)::                           i
        trace = zero
        do i=1,min(size(m,dim=1),size(m,dim=2))
            trace = trace + m(i,i)
        enddo
    end function trace
    
    
    !--bubble sort
    pure subroutine sort1_i(a,loca)
    integer(ip),dimension(:),intent(inout)::a
    integer(ip),dimension(:),&
    intent(out),optional::                  loca
    integer(ip)::                           i,k,n
    logical(lp)::                           ex
        n = size(a); ex = present(loca)
        if(ex) then
            do i=1,n
                loca(i) = i
            end do
        endif
        do i=n-1,1,-1
            do k=1,i
                if(a(k)>a(k+1)) then
                    call swap(a(k), a(k+1))
                    if(ex) call swap(loca(k), loca(k+1))
                end if
            end do
        end do
    end subroutine sort1_i
    
    !--
    pure subroutine sort2_i(a,loca)
    integer(ip),dimension(:,:),intent(inout)::  a
    integer(ip),dimension(:,:),&
    intent(out),optional::                      loca
    integer(ip)::                               j
        if(present(loca)) then
            do j=1,size(a,dim=2)
                call sort(a(:,j), loca(:,j))
            end do
        else
            do j=1,size(a,dim=2)
                call sort(a(:,j))
            end do
        endif
    end subroutine sort2_i
    
    !---
    pure subroutine sort1_r(a,loca)
    real(rp),dimension(:),intent(inout)::   a
    integer(ip),dimension(:),&
    intent(out),optional::                  loca
    integer(ip)::                           i,k,n
    logical(lp)::                           ex
        n = size(a); ex = present(loca)
        if(ex) then
            do i=1,n
                loca(i) = i
            end do
        endif
        do i=n-1,1,-1
            do k=1,i
                if(a(k)>a(k+1)) then
                    call swap(a(k), a(k+1))
                    if(ex) call swap(loca(k), loca(k+1))
                end if
            end do
        end do
    end subroutine sort1_r
    
    !--
    pure subroutine sort2_r(a,loca)
    real(rp),dimension(:,:),intent(inout):: a
    integer(ip),dimension(:,:),&
    intent(out),optional::                  loca
    integer(ip)::                           j
        if(present(loca)) then
            do j=1,size(a,dim=2)
                call sort(a(:,j), loca(:,j))
            end do
        else
            do j=1,size(a,dim=2)
                call sort(a(:,j))
            end do
        endif
    end subroutine sort2_r
    
    !-----
    pure function cumprod(x) 
    real(rp),dimension(:),intent(in)::  x
    real(rp),dimension(size(x))::       cumprod
    integer(ip)::                       i
        cumprod(1) = x(1)
        do i=2,size(x)
            cumprod(i) = cumprod(i-1)*x(i)  
        end do   
    end function cumprod
    
    !transpose matrix
    pure function trans_ivm(v) result(mt)
    integer(ip),dimension(:),intent(in)::v
    integer(ip),dimension(1,size(v))::  mt
        mt(1,:) = v
    end function trans_ivm
    !--
    pure function trans_rvm(v) result(mt)
    real(rp),dimension(:),intent(in)::  v
    real(rp),dimension(1,size(v))::     mt
        mt(1,:) = v
    end function trans_rvm
    !--
    pure function trans_rmm(m) result(mt)
    real(rp),dimension(:,:),intent(in)::    m
    real(rp),dimension(size(m,2),size(m,1))::mt
        mt = transpose(m)
    end function trans_rmm
    
    !----------------------------------------------
    pure subroutine repmat_rv(v,m,n,rep)
    real(rp),dimension(:),intent(in)::                  v
    integer(ip),intent(in)::                            m,n
    real(rp),dimension(:,:),allocatable,intent(out)::   rep
    integer(ip)::                                       i,j,di,dj
        allocate(rep(m*size(v),n))
        do j=1,n
            dj = j-1
            do i=1,m
                di = (i-1) * size(v)
                rep(di+1:di+size(v),dj+1) = v
            enddo
        enddo
    end subroutine repmat_rv
    !--
    pure subroutine repmat_rm(mat,m,n,rep)
    real(rp),dimension(:,:),intent(in)::    mat
    integer(ip),intent(in)::                m,n
    real(rp),dimension(:,:),allocatable,intent(out)::rep
    integer(ip)::                           i,j,di,dj,mm,mn
        mm = size(mat,1); mn = size(mat,2)
        allocate(rep(m*mm,n*mn))
        do j=1,n
            dj = (j-1) * mn
            do i=1,m
                di = (i-1) * mm
                rep(di+1:di+mm,dj+1:dj+mn) = mat
            enddo
        enddo
    end subroutine repmat_rm
    !--
    pure subroutine repmat_iv(v,m,n,rep)
    integer(ip),dimension(:),intent(in)::               v
    integer(ip),intent(in)::                            m,n
    integer(ip),dimension(:,:),allocatable,intent(out)::rep
    integer(ip)::                                       i,j,di,dj
        allocate(rep(m*size(v),n))
        do j=1,n
            dj = j-1
            do i=1,m
                di = (i-1) * size(v)
                rep(di+1:di+size(v),dj+1) = v
            enddo
        enddo
    end subroutine repmat_iv
    !--
    pure subroutine repmat_im(mat,m,n,rep)
    integer(ip),dimension(:,:),intent(in):: mat
    integer(ip),intent(in)::                m,n
    integer(ip),dimension(:,:),allocatable,intent(out)::rep
    integer(ip)::                           i,j,di,dj,mm,mn
        mm = size(mat,1); mn = size(mat,2)
        allocate(rep(m*mm,n*mn))
        do j=1,n
            dj = (j-1) * mn
            do i=1,m
                di = (i-1) * mm
                rep(di+1:di+mm,dj+1:dj+mn) = mat
            enddo
        enddo
    end subroutine repmat_im
    
    !--
    pure function bondcol_vvr(a,b) result(r)
    real(rp),dimension(:),intent(in)::      a,b
    real(rp),dimension(:,:),allocatable::   r
        allocate(r(max(size(a),size(b)),2))
        r = 0._rp
        r(:,1) = a
        r(:,2) = b
    end function bondcol_vvr
    
    pure function bondcol_vvvr(a,b,c) result(r)
    real(rp),dimension(:),intent(in)::      a,b,c
    real(rp),dimension(:,:),allocatable::   r
        allocate(r(max(size(a),size(b),size(c)),3))
        r = 0._rp
        r(:,1) = a
        r(:,2) = b
        r(:,2) = c
    end function bondcol_vvvr
    
    !--
    pure function bondrow_vvr(a,b) result(r)
    real(rp),dimension(:),intent(in)::      a,b
    real(rp),dimension(:,:),allocatable::   r
        allocate(r, source=trans(bondcol(a,b)))
    end function bondrow_vvr
    
    pure function bondrow_vvvr(a,b,c) result(r)
    real(rp),dimension(:),intent(in)::      a,b,c
    real(rp),dimension(:,:),allocatable::   r
        allocate(r, source=trans(bondcol(a,b,c)))
    end function bondrow_vvvr
    
    !================================================================================
    !computes the compositions of the integer n into k parts | dim(a)=k; sum(a)=|a|=n
    !e.g. n=3, k=2: (3,0)(2,1)(1,2)(0,3) | totally 4 kinds, and this subroutine accomplish this order for next
    !On the first call to this routine, set MORE = FALSE.  The routine will compute the first element in the 
    !sequence of compositions and return it, as well as setting MORE = TRUE.  If more compositions are desired
    !call again, and again.  Each time, the routine will return with a new composition.
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/subset/subset.f90 |comp_next
    pure subroutine compositionNext(n,k,a,more,h,t)
    integer(ip),intent(in)::                    n,k
    integer(ip),intent(inout)::                 h,t
    integer(ip),dimension(k),intent(inout)::    a
    logical(lp),intent(inout)::                 more
            
        if(.not.more) then
            t = n; h = 0
            a = 0; a(1) = n
        else
            !if a(1) is positive, set h=0
            if(t>1) h = 0
            !h is the entry increased last time
            !add all but one of the value to a(1) and increamenting a(h+1) by 1
            h = h + 1
            t = a(h)
            a(h) = 0
            a(1) = t - 1
            a(h+1) = a(h+1) + 1
        endif
        more = a(k) /= n
        
    end subroutine compositionNext
    
    !--
    !select k element from n element without order | no order means one order by number
    !dim(a)=k, maxval(a)<=n, no same element in a | n>k
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/subset/subset.f90 |comb_next
    pure subroutine combinationNext(n,k,a,more)
    integer(ip),intent(in)::                n,k
    integer(ip),dimension(k),intent(inout)::a
    logical(lp),intent(inout)::             more
    integer(ip)::                           i,j
        
        if(.not.more) then
            do i=1,k
                a(i) = i
            enddo
            more = n/=k
        else
            if(a(k)<n) then
                a(k) = a(k) + 1
                more = k/=1 .or. a(k)/=n
            else
                do i=k,2,-1
                    if(a(i-1)<n-k+i-1) then
                        a(i-1) = a(i-1) + 1
                        do j=i,k
                            a(j) = a(i-1) + j - (i-1)
                        enddo
                        more = a(1)/=n-k+1
                        return
                    endif
                enddo
            endif
        endif

    end subroutine combinationNext
    
    !-------------------------------------
    !    the vectors are produced in colexical order, starting with
    !    (0,        0,        ...,0),
    !    (1,        0,        ...,0),
    !     ...
    !    (base(1)-1,0,        ...,0)
    !    (0,        1,        ...,0)
    !    (1,        1,        ...,0)
    !    ...
    !    (base(1)-1,base(2)-1,...,base(dim_num)-1).
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/subset/subset.f90 %vec_colex_next2
    pure subroutine colexNext_anisotropy(k,base,a,more)
    integer(ip),intent(in)::                k
    integer(ip),dimension(:),intent(in)::   base
    integer(ip),dimension(:),intent(inout)::a
    logical(lp),intent(inout)::             more
    integer(ip)::                           i,j
        if(.not.more) then
            a = 0; more = .true.
        else
            do i=1,k
                a(i) = a(i) + 1
                if(a(1)==base(1)-1) then
                    do j=k,1,-1
                        if(a(j)/=base(j)-1) exit
                        if(j==1) more = .false.
                    enddo
                endif
                if(a(i)<base(i)) return
                a(i) = 0
            enddo
        endif
    end subroutine colexNext_anisotropy
    
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/subset/subset.f90 %vec_colex_next
    pure subroutine colexNext_isotropy(k,base,a,more)
    integer(ip),intent(in)::                k
    integer(ip),intent(in)::                base
    integer(ip),dimension(:),intent(inout)::a
    logical(lp),intent(inout)::             more
    integer(ip)::                           i,j
        if(.not.more) then
            a = 0; more = .true.
        else
            do i=1,k
                a(i) = a(i) + 1
                if(a(1)==base-1) then
                    do j=k,1,-1
                        if(a(j)/=base-1) exit
                        if(j==1) more = .false.
                    enddo
                endif
                if(a(i)<base) return
                a(i) = 0
            enddo
        endif
    end subroutine colexNext_isotropy
    
    
    !=======================================================================
    !given csr sparse matrix, give diag index of ia(j)-ia(j+1) for each row
    !ia(n+1), ja(nnz), da(n)
    !tip: ia(1) = 1, ia(n+1) = nnz + 1
    pure subroutine diagIxCsr(ia,ja,da)
    integer(ip),dimension(:),intent(in)::   ia,ja
    integer(ip),dimension(:),intent(out)::  da
    integer(ip)::                           i,j
        da = -1
        do j=1,size(da)
            do i=ia(j),ia(j+1) - 1
                if(ja(i)==j) da(j) = i
            enddo
        enddo
    end subroutine diagIxCsr
    
    
    !=======================inner product==================
    pure real(rp) function vvinnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:),intent(in)::  lhs,rhs
        p = dot_product(lhs,rhs)
    end function vvinnerproduct
    
    !--
    pure function mminnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:,:),intent(in)::lhs,rhs
    real(rp),dimension(size(lhs,dim=1),size(rhs,dim=2))::p
        p = matmul(lhs,rhs)
    end function mminnerproduct
    
    !--
    pure function mvinnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:,:),intent(in)::lhs
    real(rp),dimension(:),intent(in)::  rhs
    real(rp),dimension(size(lhs,dim=1))::p
        p = matmul(lhs,rhs)
    end function mvinnerproduct
    
    !--
    pure function vminnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:),intent(in)::  lhs
    real(rp),dimension(:,:),intent(in)::rhs
    real(rp),dimension(size(rhs,dim=2))::p
        p = matmul(lhs,rhs)
    end function vminnerproduct
    
    
    
    !---outer product
    pure function vvouterproduct(lhs,rhs) result(p)
    real(rp),dimension(:),intent(in)::  lhs,rhs
    real(rp),dimension(size(lhs),size(rhs)):: p
    integer(ip)::                       i,j
        do j=1,size(rhs)
            do i=1,size(lhs)
                p(i,j) = lhs(i) * rhs(j)
            enddo
        enddo
    end function vvouterproduct
    
    
    !---cross product
    pure function vvcrossproduct3(lhs,rhs) result(p)
    real(rp),dimension(3),intent(in)::  lhs,rhs
    real(rp),dimension(3)::             p
        p(1)  =   lhs(2) * rhs(3)   -   lhs(3) * rhs(2)
        p(2)  =   lhs(3) * rhs(1)   -   lhs(1) * rhs(3)
        p(3)  =   lhs(1) * rhs(2)   -   lhs(2) * rhs(1)
    end function vvcrossproduct3
    
    pure real(rp) function vvcrossproduct2(lhs,rhs) result(p)
    real(rp),dimension(2),intent(in)::  lhs,rhs
        p = lhs(1) * rhs(2) - lhs(2) * rhs(1)
    end function vvcrossproduct2
    
    
    !some special plus
    !--
    pure function mvplus(lhs,rhs) result(m)
    real(rp),dimension(:,:),intent(in)::lhs
    real(rp),dimension(:),intent(in)::  rhs
    real(rp),dimension(size(lhs,dim=1),size(lhs,dim=2))::m
    integer(ip)::                       i,j,d
        m = lhs
        d = size(rhs)
        forall(i=1:d,j=1:d,i==j) m(i,j) = m(i,j) + rhs(i)
    end function mvplus
    
    !--
    pure function vmplus(lhs,rhs) result(m)
    real(rp),dimension(:),intent(in)::  lhs
    real(rp),dimension(:,:),intent(in)::rhs
    real(rp),dimension(size(rhs,dim=1),size(rhs,dim=2))::m
        m = rhs + lhs
    end function vmplus
    
    !--
    pure function mvminus(lhs,rhs) result(m)
    real(rp),dimension(:,:),intent(in)::lhs
    real(rp),dimension(:),intent(in)::  rhs
    real(rp),dimension(size(lhs,dim=1),size(lhs,dim=2))::m
        m = lhs + ( - rhs )
    end function mvminus
    
    !--
    pure function vmminus(lhs,rhs) result(m)
    real(rp),dimension(:),intent(in)::  lhs
    real(rp),dimension(:,:),intent(in)::rhs
    real(rp),dimension(size(rhs,dim=1),size(rhs,dim=2))::m
        m = ( - rhs ) + lhs
    end function vmminus

    
    !equal operator
    !--
    pure logical(lp) function eqv_r(lhs,rhs) result(r)
    real(rp),dimension(:),intent(in)::      lhs,rhs
        r = .false.
        if(size(lhs)==size(rhs)) r = all(lhs==rhs)
    end function eqv_r
    
    !--
    pure logical(lp) function eqv_i(lhs,rhs) result(r)
    integer(ip),dimension(:),intent(in)::   lhs,rhs
        r = .false.
        if(size(lhs)==size(rhs)) r = all(lhs==rhs)
    end function eqv_i
    
    !--
    pure logical(lp) function eqvs(lhs,rhs) result(r)
    integer(ip),dimension(:),intent(in)::   lhs,rhs
    integer(ip),dimension(size(lhs))::      sl,sr
        r = .false.
        if(size(lhs)==size(rhs)) then
            sl = lhs; sr = rhs
            call sort(sl)
            call sort(sr)
            r = all(sl==sr)
        endif
    end function eqvs
    
end module arrayOpsLib