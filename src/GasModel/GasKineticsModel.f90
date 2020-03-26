module GasKinetics_
use constants
use arrayOpsLib
implicit none

    private
    public:: GasKinetics
    
    type GasKinetics
        
    contains
	
	
		!Sutherland Equation to calculate viscocity
        generic::               Sutherland => Sutherland_custom,Sutherland_air
        procedure,nopass::      Sutherland_custom
        procedure,nopass::      Sutherland_air
	
		!--
		procedure,nopass::		stress
		procedure,nopass::		newtonStress
        
    end type GasKinetics
        
contains
    
	!--
    elemental real(rp) function Sutherland_custom(T,T0,T1,Miu0) result(miu)
    real(rp),intent(in)::   T,T0,T1,Miu0
	
        miu = (T0 + T1)/(T + T1)
        miu = miu*Miu0*sqrt((T/T0)**3)
		
    end function Sutherland_custom
    
    !--
    elemental real(rp) function Sutherland_air(T) result(miu)
    real(rp),intent(in)::   T
	
        miu = Sutherland_custom(T, 273.15_rp, 110.4_rp, 1.7161e-5_rp)
		
    end function Sutherland_air

	!--
	function stress(miu, gradu) result(tau)
	real(rp),intent(in)::					miu
	real(rp),dimension(:,:),intent(in)::	gradu
	real(rp),dimension(size(gradu,1),&
	size(gradu,2))::						tau
	
		tau = newtonStress(miu, gradu)
	
	end function stress
	
	!--
	function newtonStress(miu, gradu) result(tau)
	real(rp),intent(in)::					miu
	real(rp),dimension(:,:),intent(in)::	gradu
	real(rp),dimension(size(gradu,1),&
	size(gradu,2))::						tau
	real(rp),dimension(size(gradu,1))::		delta
	
		delta = 1._rp
		
		tau = (gradu + transpose(gradu))/2._rp
        tau = 2._rp*miu*(tau - 1._rp/3._rp*trace(gradu)*delta)
	
	end function newtonStress

end module GasKinetics_