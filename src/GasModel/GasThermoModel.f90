module GasThermo_
use constants
implicit none

    private
    public:: GasThermo
    
    type GasThermo
        
    contains
    
		procedure,nopass::			FourierLaw
	
        !-------------
        !idealGasModel
        !-------------
        !ideal Gas state Equation T(p,rho,R)
        procedure,nopass::          T_ideal

        !ideal gas isoEntropy flow
        procedure,nopass::          Pr_s => isoEntropyPressureRatio
        procedure,nopass::          Tr_s => isoEntropyTemperatureRatio
        
    end type GasThermo
        
contains
    
    function FourierLaw(k, gradT) result(q)
	real(rp),intent(in)::				k
	real(rp),dimension(:),intent(in)::	gradT
	real(rp),dimension(size(gradT))::	q
	
		q = -k*gradT
	
 	end function FourierLaw
    
	!--
    elemental real(rp) function T_ideal(rho,p,R) result(T)
    real(rp),intent(in)::   rho,p,R
	
        T = p/rho/R
		
    end function T_ideal
    
    !p0/p for ideal isoEntropy flow
    elemental real(rp) function isoEntropyPressureRatio(ma,gm) result(pr)
    real(rp),intent(in)::   ma,gm
	
        pr = isoEntropyTemperatureRatio(ma,gm)**(gm/(gm - 1._rp))
		
    end function isoEntropyPressureRatio
    
    !T0/T for ideal isoEntropy flow
    elemental real(rp) function isoEntropyTemperatureRatio(ma,gm) result(TR)
    real(rp),intent(in)::   ma,gm
	
        TR = 1._rp + 0.5_rp*(gm - 1._rp)*ma**2
		
    end function isoEntropyTemperatureRatio

end module GasThermo_