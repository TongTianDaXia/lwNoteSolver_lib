module sensitivity_
use constants
use stringOpsLib
use IntegrationLib
use SpecialFunctionLib
use polynomial_
use multiIndex_
implicit none

    private
    public::    sensitivity
    
    !--
    type:: sensitivity
    
        private
        type(polynomial),dimension(:),pointer:: basis_
        type(multiAnovaIndex)::                 ind_
        integer(ip)::                           p_
        integer(ip)::                           dim_
        
    contains
    
        generic::       init => init_0,init_basis
        procedure::     init_0
        procedure::     init_basis
        
        !Global sensitivity based on deviation
        procedure::     Global_dev
        
    end type
    
contains
    
    pure subroutine init_0(this)
    class(sensitivity),intent(out)::  this
    end subroutine init_0
    
    !p the order of polynomial | p = ubound(basis)
    !d the dimenson
    subroutine init_basis(this,basis,p,d)
    class(sensitivity),intent(out)::  this
    type(polynomial),dimension(0:),target,intent(in)::  basis
    integer(ip),intent(in)::                            p,d
        this%basis_ => basis
        this%p_ = p
        this%dim_ = d
        call this%ind_%init(p, d)
    end subroutine init_basis
    
    !main subroutine | computing sensitivity of different combination of dimension
    subroutine Global_dev(this,f,quadx,quadw,QuadratureRule,senIndt,senInd1,senInd2,senInd3)
    class(sensitivity),intent(in)::                     this
    real(rp),dimension(:),intent(in)::                  f               !f(x at np), f value at i-th point of x
    real(rp),dimension(:,:),intent(in)::                quadx           !quadrature x
    real(rp),dimension(:),intent(in)::                  quadw           !quadrature w
    character(*),intent(in)::                           QuadratureRule  !quadrature rule
    real(rp),dimension(:),intent(out)::                 senIndt         !total sensitivity index
    real(rp),dimension(:),intent(out)::                 senInd1         !first order sensitivity index
    real(rp),dimension(:,:),intent(out)::               senInd2         !second order sensitivity index
    real(rp),dimension(:,:,:),intent(out),optional::    senInd3         !third order sensitivity index
    integer(ip)::                                       i,j,k,np,indtn              
    integer(ip)::                                       dimSg,npSg      !sg dimension and sg sample number
    integer(ip),dimension(:),allocatable::              dimr                        
    integer(ip),dimension(:,:),allocatable::            alphaInd        !alpha set
    real(rp)::                                          scale,m,dev                 
    real(rp),dimension(:),allocatable::                 cubaWeight_
    character(len(trim(quadratureRule)))::              rule
    
        !----------------------------
        npSg = size(f, 1)
        dimSg = this%dim_
        np = merge(3, 2, present(senInd3))
        
        allocate(dimr(np))
        do i=1,np                                       
            dimr(i) = this%ind_%idxdim(i)
        end do 
        indtn = sum(dimr(:))
        
        !--
        allocate(alphaInd(indtn, dimSg))
        allocate(cubaWeight_(indtn))
        alphaInd = 0
        cubaWeight_ = 0._rp
        
        !--
        !quadratureRule = lowerString(quadratureRule)
        rule = trim(quadratureRule)
        call lowerString(rule)
        if(rule=='gausshermite' .or. rule=='gh') then
            scale = srpi
        elseif(rule=='clenshawcurtis' .or. rule=='cc') then
            scale = 2._rp
        elseif(rule=='gausslegendre' .or. rule=='gl') then
            scale = 2._rp
        endif
        
        !----------------------------
        m = sum(quadw(1:npSg)*f(1:npSg))/(scale**dimSg)
        dev = sum(quadw(1:npSg)*(f(1:npSg) - m)**2)/(scale**dimSg)
        if(dev<0._rp) then
            print*,'error: Sensitivity/Global_dev has unacceptatble deviation with value less than 0'
            stop
        end if
        
        !----------------------------
        call alphaSet(alphaInd)
        call cubeWeight(cubaWeight_)
        call senIndex1st(senInd1)
        call senIndex2nd(senInd2)
        
        !--
        if(present(senInd3)) then 
            call senIndex3rd(senInd3)
            call senIndexTot(senIndt, senInd1, senInd2, senInd3)
        else
            call senIndexTot(senIndt, senInd1, senInd2)
        end if
        
    contains
        
        pure subroutine alphaSet(alphaI)
        integer(ip),dimension(:,:),intent(out)::    alphaI
        integer(ip)::                               i,j,k,alpn
        logical(lp)::                               more

            alpn = 1
            do i=1,np
                more = .true.
                do j=alpn, alpn + dimr(i) - 1
                    call this%ind_%traverse(i,alphaI(j,:), more)
                end do
                alpn = sum(dimr(1:i)) + 1 
            end do 
            
        end subroutine alphaSet
        
        !-----------------------------
        pure subroutine cubeWeight(cw)
        real(rp),dimension(:),intent(out)::         cw
        integer(ip)::                               i,j,k,alpn
        real(rp),dimension(:),allocatable::         mPolyValue
        
            allocate(mPolyValue(npSg))
            alpn = 1
            do i=1,np
                do j=alpn,alpn + dimr(i) - 1
                    do k=1,npSg
                        mPolyValue(k) = multiPolynominal(this%basis_, alphaInd(j,:), quadx(:,k))
                    end do
                    cw(j) = sum(quadw(:)*f(1:npSg)*mPolyValue(:))/(scale**dimSg)
                end do
                alpn = sum(dimr(1:i)) + 1
            end do
            
        end subroutine cubeWeight
        
        !-----------------------------
        subroutine senIndex1st(s1)
        real(rp),dimension(:),intent(out)::         s1
        integer(ip)::                               i,j,k,z
        real(rp),dimension(:),allocatable::         mPolyValue
        real(rp)::                                  mpolyDev
            s1(:) = 0._rp
            allocate(mPolyValue(npSg))
            do z=1,dimSg
                do j=1,dimr(1)
                    if(alphaInd(j,z)/=0) then
                        do k=1,npSg
                            mPolyValue(k) = multiPolynominal(this%basis_, alphaInd(j,:), quadx(:,k))
                        end do
                        mpolyDev = sum(quadw(:)*mPolyValue(:)*mPolyValue(:))/(scale**dimSg)
                        s1(z) = s1(z) + mpolyDev*cubaWeight_(j)*cubaWeight_(j)
                    endif
                enddo
            enddo
            s1(:) = s1(:)/dev
        end subroutine senIndex1st
        
        !-----------------------------
        pure subroutine senIndex2nd(s2)
        real(rp),dimension(:,:),intent(out)::       s2
        integer(ip)::                               i,j,k,z1,z2
        real(rp),dimension(:),allocatable::         mPolyValue
        real(rp)::                                  mpolyDev
            s2(:,:) = 0._rp 
            allocate(mPolyValue(npSg))
            do z1=1,dimSg-1
                do z2=z1+1,dimSg
                    do j=dimr(1)+1,dimr(1)+dimr(2)
                        if(alphaInd(j,z1)/=0.and.alphaInd(j,z2)/=0) then
                            do k=1,npSg
                                mPolyValue(k) = multiPolynominal(this%basis_, alphaInd(j,:), quadx(:,k))  
                            end do
                            mpolyDev = sum(quadw(:)*mPolyValue(:)*mPolyValue(:))/(scale**dimSg)
                            s2(z1,z2) = s2(z1,z2) + mpolyDev*cubaWeight_(j)*cubaWeight_(j)
                        endif
                    enddo
                enddo
            enddo
            s2(:,:) = s2(:,:)/dev
        end subroutine senIndex2nd
        
        !-----------------------------
        pure subroutine senIndex3rd(s3)
        real(rp),dimension(:,:,:),intent(out)::     s3
        integer(ip)::                               i,j,k,z1,z2,z3
        real(rp),dimension(:),allocatable::         mPolyValue
        real(rp)::                                  mpolyDev
            s3(:,:,:) = 0._rp
            allocate(mPolyValue(npSg))
            do z1=1,dimSg-2
                do z2=z1+1,dimSg-1
                    do z3=z2+1,dimSg
                        do j=dimr(1)+dimr(2)+1,dimr(1)+dimr(2)+dimr(3)
                            if(alphaInd(j,z1)/=0.and.alphaInd(j,z2)/=0.and.alphaInd(j,z3)/=0) then
                                do k=1,npSg
                                    mPolyValue(k) = multiPolynominal(this%basis_, alphaInd(j,:), quadx(:,k))
                                end do
                                mpolyDev = sum(quadw(:)*mPolyValue(:)*mPolyValue(:))/(scale**dimSg)
                                s3(z1,z2,z3) = s3(z1,z2,z3) + mpolyDev*cubaWeight_(j)*cubaWeight_(j)
                            end if
                        end do
                    end do
                end do
            end do
            s3(:,:,:) = s3(:,:,:)/dev
        end subroutine senIndex3rd
        
        !--
        pure subroutine senIndexTot(st,s1,s2,s3)
        real(rp),dimension(:),intent(out)::             st
        real(rp),dimension(:),intent(in)::              s1
        real(rp),dimension(:,:),intent(in)::            s2
        real(rp),dimension(:,:,:),intent(in),optional:: s3
        integer(ip)::                                   i,j,k,z
            
            st(:) = s1(:)
            !--
            do i=1,dimSg
                do j=1,dimSg-1
                    do k=j+1,dimSg
                        if(j==i.or.k==i) st(i) = st(i) + s2(j,k)
                    enddo
                enddo
            enddo
            
            !--
            if(present(s3)) then
                do i=1,dimSg
                    do j=1,dimSg-2
                        do k=j+1,dimSg-1
                            do z=k+1,dimSg
                                if(j==i.or.k==i.or.z==i) st(i) = st(i) + s3(j,k,z)
                            enddo
                        enddo
                    enddo
                enddo
            endif
        
        end subroutine senIndexTot 
        
    end subroutine Global_dev
    
end module sensitivity_