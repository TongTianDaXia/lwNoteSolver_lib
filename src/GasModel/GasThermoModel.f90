module GasThermo_
use constants
implicit none

    private
    public:: GasThermo
    
    type GasThermo
        
    contains
    
        !Sutherland Equation to calculate viscocity
        generic::                   Sutherland => Sutherland_custom,Sutherland_air
        procedure,nopass::          Sutherland_custom
        procedure,nopass::          Sutherland_air
        
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
    
    !---------------
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
    
!------------------------------------------------------------------------
    elemental real(rp) function T_ideal(rho,p,R) result(T)
    real(rp),intent(in)::   rho,p,R
        T = p / rho / R
    end function T_ideal
    
    !p0/p for ideal isoEntropy flow
    elemental real(rp) function isoEntropyPressureRatio(ma,gm) result(pr)
    real(rp),intent(in)::   ma,gm
        pr = isoEntropyTemperatureRatio(ma,gm)**(gm/(gm-1._rp))
    end function isoEntropyPressureRatio
    
    !T0/T for ideal isoEntropy flow
    elemental real(rp) function isoEntropyTemperatureRatio(ma,gm) result(TR)
    real(rp),intent(in)::   ma,gm
        TR = 1._rp + 0.5_rp*(gm-1._rp)*ma**2
    end function isoEntropyTemperatureRatio

end module GasThermo_
