module GasSpecies_
use constants
use stringOpsLib
implicit none

	private
	public:: GasSpecies

	!--
    type GasSpecies

        private

		character(16)::			name_

		!criterion temperature for thermodynamic fit function
		real(rp),dimension(:),&
		allocatable::			Tc_
        real(rp),dimension(:,:),&
        allocatable::           thermoCoef_ 	!(ncoef, ix) | [ix = 1:size(Tc_)-1]

		!--
        real(rp)::              mw_         	![g/mol]
        real(rp)::              R_          	!J/[K*g] = R_c / mw

		!--
		!Lennard-Jones potential: V(r) = 4\epsilon((\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6)
		integer(ip)::			geo_			!geometry 
        real(rp)::              T0_   			!T0 = epsilon/k [K]
												!k: boltzmann constant
                                                !epsilon: Lennard-Jones potential well depth
        real(rp)::              sigma_          !Lennard-Jones collision diameter [Angstroms] = 1.e-10[meter]
		real(rp)::				muDip_			!dipole moment [Debye] = 1.e-18[cm^{3/2}erg^{1/2}]
		!note:: the unit of Debye in chemkin seems different from the ohter definition
		!~ generally [Debye] =  3.33564e¨C30 [C m] = 3.33564e¨C30 [A s m]
		real(rp)::				alpha_			!polarizabilities [Angstroms] = 1.e-10[meter]
		real(rp)::				Zrot298_		!rotational relaxation collision number

    contains

        generic::				init    =>  init_n, &
								            init_spinfo, &
											init_spinfo_2T
        procedure::				init_n
        procedure::				init_spinfo
		procedure::				init_spinfo_2T

        !menber function
        generic::				thermoCoef => thermoCoef_ij, &
												thermoCoef_ptrj, &
												thermoCoef_T
        procedure::				thermoCoef_ij
        procedure::				thermoCoef_ptrj
		procedure::				thermoCoef_T

		!--
        procedure::				mw => molecularWeight
        procedure::				MolecularWeight ![g/mol]

        procedure::				R               ![J/(K*g)] | R_c/mw

        !based on the
        generic::				Cp => CpNasa9_R,CpNasa9_Rc
        procedure::				CpNasa9_R	!dimension same as R, if R[J/(K*g)], then output [J/(K*g)]
        procedure::				CpNasa9_Rc	![J/(K*mol)], R_c as the default R

        !H = \int_{0}^{T} Cp d{T}
        generic::               H => EnthalpyNasa9_R, EnthalpyNasa9_Rc
        procedure::             EnthalpyNasa9_R		!dimension same as R, if R[J/(K*g)], then output [J/g]
        procedure::             EnthalpyNasa9_Rc	![J/mol], R_c as the default R

        !S = \int_{0}^{T} \frac{Cp d{T}}{T}
        generic::               S => EntropyNasa9_R, EntropyNasa9_Rc
        procedure::             EntropyNasa9_R		!dimension same as R, if R[J/(K*g)], then output [J/g]
        procedure::             EntropyNasa9_Rc		![J/(K*mol)], R_c as the default R

        !transport parameter
        procedure::             mu => viscosity_sp	![kg/m/s]
		procedure::				D => diffusion_sp	![m^2/s]
		
		!menber function
		procedure::				geo
		procedure::				T0	
		procedure::				sigma
		procedure::				muDip
		procedure::				alpha
    	procedure::				Zrot298
		
	end type GasSpecies
	

contains

	!-------------
	pure subroutine init_n(this, n, m)
	class(GasSpecies),intent(out)::    	this
	integer(ip),intent(in)::            n, m

		this%name_ = 'none'
	    allocate(this%thermoCoef_(n, m))
	    this%thermoCoef_ = 0._rp
	    this%mw_ = 0._rp
	    this%R_ = 0._rp
	    this%T0_ = 0._rp
	    this%sigma_ = 0._rp

	end subroutine init_n

	!--mw[g/mol]
	subroutine init_spinfo(this, name, mw, Tc, thermoCoef, T0, sigma)
	class(GasSpecies),intent(out)::			this
	character(*),intent(in)::				name
	real(rp),intent(in)::               	mw
	real(rp),dimension(:,:),optional,&
	intent(in)::							thermoCoef
	real(rp),dimension(:),optional,&
	intent(in)::							Tc
	real(rp),optional,intent(in)::			T0, sigma

		if(present(Tc).and.present(thermoCoef)) then
			if(size(Tc)/=size(thermoCoef, 2) + 1) &
			stop 'error: GasSpeciesModel/init_spinfo get a dismatching dimension of Tc and thermoCoef'
		endif

		this%name_ = lowerStringFc(name)
		this%mw_ = mw
	    this%R_ = R_c/mw
		
		if(present(Tc)) allocate(this%Tc_, source=Tc)
	    if(present(thermoCoef)) allocate(this%thermoCoef_, source=thermoCoef)
	    
		if(present(T0)) then
		    this%T0_ = T0
		else
			this%T0_ = 0._rp
		endif

		if(present(sigma)) then
		    this%sigma_ = sigma
		else
			this%sigma_ = 0._rp
		endif

	end subroutine init_spinfo

	!--mw[g/mol]
	subroutine init_spinfo_2T(this, name, mw, Tc, thermoCoef, T0, sigma)
	class(GasSpecies),intent(out)::			this
	character(*),intent(in)::				name
	real(rp),intent(in)::               	Tc, mw
	real(rp),dimension(:,:),intent(in)::	thermoCoef
	real(rp),optional,intent(in)::			T0, sigma

		if(present(T0) .and. present(sigma)) then
			call this%init(name, mw, [250._rp, Tc, 5000._rp], thermoCoef, T0, sigma)
		else
			call this%init(name, mw, [250._rp, Tc, 5000._rp], thermoCoef)
		endif

	end subroutine init_spinfo_2T

	!--------------------
	elemental real(rp) function thermoCoef_ij(this, i, j) result(c)
	class(GasSpecies),intent(in)::          this
	integer(ip),intent(in)::                i, j

	    c = this%thermoCoef_(i, j)

	end function thermoCoef_ij

	!--
	function thermoCoef_ptrj(this, j) result(p)
	class(GasSpecies),target,intent(in)::   this
	integer(ip),intent(in)::				j
	real(rp),pointer,dimension(:)::         p

	    p => this%thermoCoef_(:, j)

	end function thermoCoef_ptrj

	!--
	function thermoCoef_T(this, T) result(p)
	class(GasSpecies),target,intent(in)::   this
	real(rp),intent(in)::					T
	real(rp),pointer,dimension(:)::         p
	integer(ip)::							i

		if(T<this%Tc_(1)) then
			stop 'error: GasChemModel/thermoCoef_T get a too low temperature'
		else
			do i=2,size(this%Tc_)
				if(T<this%Tc_(i)) exit
			enddo
		endif

		if(i>size(this%Tc_)) stop 'error: GasChemModel/thermoCoef_T get a too high temperature'

		p => this%thermoCoef(i-1)

	end function thermoCoef_T

	!-------
	elemental real(rp) function molecularWeight(this)
	class(GasSpecies),intent(in)::     this

	    MolecularWeight = this%mw_

	end function molecularWeight

	!-------
	elemental real(rp) function R(this)
	class(GasSpecies),intent(in)::     this

	    R = this%R_

	end function R

	!specific heat capacity polynomial fitting
	impure elemental real(rp) function CpNasa9_R(this, T, R) result(cp)
	class(GasSpecies),intent(in)::      this
	real(rp),intent(in)::               T, R
	real(rp),dimension(:),pointer::     c

		c => this%thermoCoef(T)
	    cp = R*sum(c(1:7)*T**[-2:4])

	end function CpNasa9_R

	impure elemental real(rp) function CpNasa9_Rc(this, T) result(cp)
	class(GasSpecies),intent(in)::      this
	real(rp),intent(in)::               T

	    cp = this%Cp(T, R_c)

	end function CpNasa9_Rc

	!-----
	!Specific heat enthalpy polynomial fitting
	impure elemental real(rp) function EnthalpyNasa9_R(this, T, R) result(enthalpy)
	class(GasSpecies),intent(in)::      this
	real(rp),intent(in)::               T, R
	real(rp),pointer,dimension(:)::     c

	    c => this%thermoCoef(T)
	    enthalpy = -c(1)/T + c(2)*log(T) + c(8)
	    enthalpy = enthalpy + sum(c(3:7)*T**[1:5]/real([1:5], rp))
	    enthalpy = R*enthalpy

	end function EnthalpyNasa9_R

	!Specific heat enthalpy polynomial fitting
	impure elemental real(rp) function EnthalpyNasa9_Rc(this, T) result(enthalpy)
	class(GasSpecies),intent(in)::      this
	real(rp),intent(in)::               T

	    enthalpy = this%H(T, R_c)

	end function EnthalpyNasa9_Rc

	!----
	!specific entropy polynomial fitting
	impure elemental real(rp) function EntropyNasa9_R(this, T, R) result(entropy)
	class(GasSpecies),intent(in):: 		this
	real(rp),intent(in)::               T, R
	real(rp),pointer,dimension(:)::     c

	    c => this%thermoCoef(T)
	    entropy = -0.5_rp*c(1)/T**2 - c(2)/T + c(3)*log(T) + c(9)
	    entropy = entropy + sum(c(4:7)*T**[1:4]/real([1:4], rp))
	    entropy = R*entropy

	end function EntropyNasa9_R

	!specific entropy polynomial fitting
	impure elemental real(rp) function EntropyNasa9_Rc(this, T) result(entropy)
	class(GasSpecies),intent(in)::  this
	real(rp),intent(in)::           T

	    entropy = this%S(T, R_c)

	end function EntropyNasa9_Rc

	!--refer to chemkin manual [pure species vicosity]
	elemental real(rp) function viscosity_sp(this, T) result(mu)
	class(GasSpecies),intent(in)::  this
	real(rp),intent(in)::           T
	real(rp)::                      omega22, Ts

	    !refer to
		![a fortran computer code package for the evaluation of gas-phase, multicomponent transport properties]
	    Ts = T/this%T0_

		!model 1 is verified with the below
	    omega22 = 1.147_rp*Ts**(-0.145_rp) + (Ts + 0.5_rp)**(-2)
		!model 2, slower than above
		!omega22 = 1.16145_rp/Ts**0.14874_rp &
		!		+ 0.52487_rp/exp(0.7732_rp*Ts) &
		!		+ 2.16178_rp/exp(2.43787_rp*Ts)

	    !mu = (5/16)[sqrt(pi*m*k*T)/pi/sigma**2/omega]
	    mu = 2.6693e-6_rp*(sqrt(this%mw_*T)/this%sigma_**2/omega22)

	end function viscosity_sp
	
	!--
	real(rp) function diffusion_sp(this, that, T, p, mod) result(D)
	class(GasSpecies),intent(in)::		this, that
	real(rp),intent(in)::				T, p
	integer(ip),optional,intent(in)::	mod
	real(rp)::							omega11, Ts, Ds, mw, T0, sigma, xi
	real(rp)::							sigmam, T0m, unitFactor
	real(rp),external::					omeg12, D12
	
		unitFactor = (1.e-18_rp)**2*(1.e8_rp)**3/(k_b*1.e7_rp)
		sigmam = 0.5_rp*(this%sigma_ + that%sigma_)
		T0m = sqrt(this%T0_*that%T0_)
		
		if(this%muDip_>tinrp .and. that%muDip_<tinrp) then
			!<this> is polar molecule
			xi = xifc(this, that)
			Ds = 0._rp
		elseif(this%muDip_<tinrp .and. that%muDip_>tinrp) then
			!<that> is polar molecule
			xi = xifc(that, this)
			Ds = 0._rp
		else
			xi = 1._rp
			Ds = 0.5_rp*this%muDip_*that%muDip_/(T0m*sigmam**3)*unitFactor
		endif
		
		mw = this%mw_*that%mw_/(this%mw_ + that%mw_)
		sigma = sigmam*xi**(-1._rp/6._rp)
		T0 = T0m*xi**2
		
		Ts = T/T0m
		D = 0.0188346_rp*sqrt(T**3/mw)/p
		
		if(present(mod)) then
			!L-J 12-6 model
			selectcase(mod)
			case(0)
				D = D/(sigmam**2*omeg12(1, Ts, Ds))
				D = D*D12(this%mw_, that%mw_, T, this%T0_, that%T0_, this%sigma_, that%sigma_, Ds)
			case(1)
				omega11 = 1.06036_rp/Ts**0.15610_rp &
						+ 0.19300_rp/exp(0.47635_rp*Ts) &
						+ 1.03587_rp/exp(1.52996_rp*Ts) &
						+ 1.76474_rp/exp(3.89411_rp*Ts)
				D = D/(sigmam**2*omega11)
			case default
				stop 'error: GasSpeciesModel%diffusion_sp get an unknown model for omega'
			endselect
		else
			D = D/(sigmam**2*omeg12(1, Ts, Ds))
			D = D*D12(this%mw_, that%mw_, T, this%T0_, that%T0_, this%sigma_, that%sigma_, Ds)
		endif
		
	contains
	
		!refer to subroutine (chemkin%tranfit.f%diffit), different from the manual of chemkin
		!eg. code write muDisp^2 and manual write muDisp for xifc
		pure real(rp) function xifc(p, np)
		class(GasSpecies),intent(in)::	p, np
		real(rp)::						alphas, muDips2, pEpsilon, pMuDip, pSigma
		
			!polar
			pEpsilon = p%T0_*k_b*1.e7_rp	![J] -> [erg]
			pMuDip = p%muDip_*1.e-18_rp		![Debye] -> [cm^{3/2}erg^{1/2}] according to Chemkin definition
			pSigma = p%sigma_*1.e-8			![angstrom] -> [cm]
			muDips2 = pMuDip**2/(pEpsilon*pSigma**3)
			
			!non-polar
			alphas = np%alpha_/np%sigma_**3 !no unit
			
			xifc = 1._rp + 0.25_rp*alphas*muDips2*sqrt(p%T0_/np%T0_)
			
		end function xifc
		
	end function diffusion_sp
	
	!--
	elemental integer(ip) function geo(this)
	class(GasSpecies),intent(in)::	this
	
		geo = this%geo_
	
	end function geo
	
	!--
	elemental real(rp) function T0(this)
	class(GasSpecies),intent(in)::	this
	
		T0 = this%T0_
	
	end function T0
	
	!--
	elemental real(rp) function sigma(this)
	class(GasSpecies),intent(in)::	this
	
		sigma = this%sigma_
	
	end function sigma
	
	!--
	elemental real(rp) function muDip(this)
	class(GasSpecies),intent(in)::	this
	
		muDip = this%muDip_
	
	end function muDip
	
	!--
	elemental real(rp) function alpha(this)
	class(GasSpecies),intent(in)::	this
	
		alpha = this%alpha_
	
	end function alpha
	
	!--
	elemental real(rp) function Zrot298(this)
	class(GasSpecies),intent(in)::	this
	
		Zrot298 = this%Zrot298_
	
	end function Zrot298


end module GasSpecies_
