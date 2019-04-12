module GasChemReactKin_
use constants, only: ip, rp, R_c, P_atm, GlobalEps
implicit none

    private          
    public:: GasChemReactKin
    public:: GasSpecies
    !--
    public:: miu_mix
    
    !----------------------------------------------------------------
    type GasSpecies
        
        private
        real(rp),dimension(:,:),&
        allocatable::               ThermoCoef_ !(ncoef,ix)|[ix=1 for low][ix=2 for high]
        real(rp)::                  Tc_         !criterion temperature for thermodynamic fit function
        real(rp)::                  mw_         ![g/mol]
        real(rp)::                  R_          !J/[K*g] = R_c / mw
        
        real(rp)::                  epsilonOverk_   !epsilon/k | k: boltzmann constant
                                                    !epsilon: Lennard-Jones potential well depth
        real(rp)::                  sigma_          !Lennard-Jones collision diameter
        
    contains  
          
        generic::                   init    =>  init_n,     &
                                                init_spinfo
        procedure::                 init_n
        procedure::                 init_spinfo

        !menber function
        generic::                   ThermoCoefH => ThermoCoefH_i,ThermoCoefH_ptr
        procedure::                 ThermoCoefH_i
        procedure::                 ThermoCoefH_ptr
        
        generic::                   ThermoCoefL => ThermoCoefL_i,ThermoCoefL_ptr
        procedure::                 ThermoCoefL_i
        procedure::                 ThermoCoefL_ptr
        
        procedure::                 mw => molecularWeight
        procedure::                 MolecularWeight ![g/mol]
        
        procedure::                 R               ![J/(K*g)] | R_c/mw
        
        !based on the 
        generic::                   Cp => CpNasa9_R,CpNasa9_Rc
        procedure::                 CpNasa9_R         !dimension same as R, if R[J/(K*g)], then output [J/(K*g)]
        procedure::                 CpNasa9_Rc        ![J/(K*mol)], R_c as the default R
        
        !H = \int_{0}^{T} Cp d{T}
        generic::                   H => EnthalpyNasa9_R, EnthalpyNasa9_Rc
        procedure::                 EnthalpyNasa9_R   !dimension same as R, if R[J/(K*g)], then output [J/g]
        procedure::                 EnthalpyNasa9_Rc  ![J/mol], R_c as the default R
        
        !S = \int_{0}^{T} \frac{Cp d{T}}{T}
        generic::                   S => EntropyNasa9_R, EntropyNasa9_Rc
        procedure::                 EntropyNasa9_R    !dimension same as R, if R[J/(K*g)], then output [J/g]
        procedure::                 EntropyNasa9_Rc   ![J/(K*mol)], R_c as the default R
        
        !miu[]
        procedure::                 miu => viscosity_sp
        
    end type GasSpecies
    
    
    !------------------------------!
    type GasChemReactKin
        
        private
        integer(ip)::                           ns_=0, nr_=0
        
        type(GasSpecies),dimension(:),&
        allocatable::                           species_

        integer(ip),dimension(:,:,:),&
        allocatable::                           stoichiometric_ !(nr,ns,ix) [ix=1 forward][ix=2 backward]
        real(rp),dimension(:,:),allocatable::   threebody_  !(nr,ns)coef for ThreebodyCon
        real(rp),dimension(:,:),allocatable::   arns_       !(3,nr) coef for arrhenius
        !--
        real(rp),dimension(:,:),allocatable::   arnsRev_!(3,nr) coef for Rev rateconstant
    
    contains
    
        generic::                   init => init_sp, init_spReact
        procedure::                 init_sp
        procedure::                 init_spReact
        procedure::                 init_movealloc
        
        !Chemistry production of basic model, see ChemKin of Ansys
        !mol generation rate[mol/cm3/s]
        generic::                   molPr => ProdRate_Basic, ProdRate_mc
        procedure,nopass::          ProdRate_Basic
        procedure::                 ProdRate_mc
        
        !mol/cm3/s
        procedure,nopass::          ProgRate
        procedure,nopass::          ProgRate1
        
        !calculate rate constant => k = a T^b \exp{- \frac{Ea}{T}}
        generic::                   k => Arrhenius
        procedure,nopass::          Arrhenius
        
        procedure::                 kf => Arrheniusf_i
        procedure::                 Arrheniusf_i
        
        !calculate Rev Rate Constant based on Equilibrium assumption
        generic::                   Revk => RevRateConstant_equilibrium
        procedure,nopass::          RevRateConstant_equilibrium
        !procedure,nopass::          RevRateConstant_Saha
        
        !calculate Equilibrium Constant based on experience formula
        procedure,nopass::          EquilibriumConstant
        
        !--
        generic::                   TbCon => ThreebodyCon,ThreebodyCon_molcon
        procedure,nopass::          ThreebodyCon
        procedure::                 ThreebodyCon_molcon
        
        !--
        generic::                   T => temperature_sp,temperature_kin
        procedure,nopass::          Temperature_sp
        procedure::                 temperature_kin
        
        !--
        procedure::                 Cp => Cp_mf     !j/(K*g)
        procedure::                 H => enthalpy_mf!j/g
        procedure::                 S => entropy_mf !j/g
        procedure::                 R => R_mf       !j/(K*g)
        
        !menber function
        generic::                   sp => sp_ptr,sp_i
        procedure::                 sp_i
        procedure::                 sp_ptr
        
        procedure::                 smf
        procedure::                 smr
        
        procedure::                 ns
        procedure::                 nr
        
    end type GasChemReactKin
    
    !-------------
        
contains

    !-------------
    pure subroutine init_n(this,n,m)
    class(GasSpecies),intent(out)::    this
    integer(ip),intent(in)::            n,m
        allocate(this%ThermoCoef_(n,m))
        this%ThermoCoef_ = 0._rp
        this%mw_ = 0._rp
        this%R_ = 0._rp
        this%epsilonOverk_ = 0._rp
        this%sigma_ = 0._rp
    end subroutine init_n
    !--mw[g/mol]
    pure subroutine init_spinfo(this, ThermoCoef, Tc, mw, epsilonOverk, sigma)
    class(GasSpecies),intent(out)::    this
    real(rp),dimension(:,:),intent(in)::ThermoCoef
    real(rp),intent(in)::               Tc, mw, epsilonOverk, sigma
        allocate(this%ThermoCoef_ , source = ThermoCoef)
        this%Tc_ = Tc
        this%mw_ = mw
        this%R_ = R_c/mw
        this%epsilonOverk_ = epsilonOverk
        this%sigma_ = sigma
    end subroutine init_spinfo
    

    !--------------------!
    elemental real(rp) function ThermoCoefL_i(this,i) result(c)
    class(GasSpecies),intent(in)::          this
    integer(ip),intent(in)::                i
        c = this%ThermoCoef_(i,1)
    end function ThermoCoefL_i
    
    function ThermoCoefL_ptr(this) result(p)
    class(GasSpecies),target,intent(in)::  this
    real(rp),pointer,dimension(:)::             p
        p => this%ThermoCoef_(:,1)
    end function ThermoCoefL_ptr
    
    !--------------------!
    elemental real(rp) function ThermoCoefH_i(this,i) result(c)
    class(GasSpecies),intent(in)::          this
    integer(ip),intent(in)::                i
        c = this%ThermoCoef_(i,2)
    end function ThermoCoefH_i
    
    function ThermoCoefH_ptr(this) result(p)
    class(GasSpecies),target,intent(in)::   this
    real(rp),pointer,dimension(:)::         p
        p => this%ThermoCoef_(:,2)
    end function ThermoCoefH_ptr
    
    
    !--------------------
    elemental real(rp) function molecularWeight(this)
    class(GasSpecies),intent(in)::     this
        MolecularWeight = this%mw_
    end function molecularWeight
    
    
    !--------------------
    elemental real(rp) function R(this)
    class(GasSpecies),intent(in)::     this
        R = this%R_
    end function R
    
    
    !specific heat capacity polynomial fitting
    elemental real(rp) function CpNasa9_R(this,T,R) result(cp)
    class(GasSpecies),intent(in)::      this
    real(rp),intent(in)::               T,R
    real(rp),dimension(:),pointer::     c
        
        if(T>this%Tc_) then
            c => this%ThermoCoefH()
        else
            c => this%ThermoCoefL()
        end if
        cp = R * sum(c(1:7)*T**[-2:4])
        
    end function CpNasa9_R
    
    elemental real(rp) function CpNasa9_Rc(this,T) result(cp)
    class(GasSpecies),intent(in)::      this
    real(rp),intent(in)::               T
        cp = this%Cp(T,R_c)
    end function CpNasa9_Rc
    
!------------------------------------------------------------------------
    !Specific heat enthalpy polynomial fitting
    elemental real(rp) function EnthalpyNasa9_R(this,T,R) result(enthalpy)
    class(GasSpecies),intent(in)::      this
    real(rp),intent(in)::               T,R
    real(rp),pointer,dimension(:)::     c
    
        if(T>this%Tc_) then
            c => this%ThermoCoefH()
        else
            c => this%ThermoCoefL()
        end if
        enthalpy = -c(1)/T + c(2)*log(T) + c(8)
        enthalpy = enthalpy + sum(c(3:7)*T**[1:5]/real([1:5],kind=rp))
        enthalpy = R * enthalpy
        
    end function EnthalpyNasa9_R
    
    !Specific heat enthalpy polynomial fitting
    elemental real(rp) function EnthalpyNasa9_Rc(this,T) result(enthalpy)
    class(GasSpecies),intent(in)::      this
    real(rp),intent(in)::               T
        enthalpy = this%H(T,R_c)
    end function EnthalpyNasa9_Rc
    
    
!------------------------------------------------------------------------
    !specific entropy polynomial fitting
    elemental real(rp) function EntropyNasa9_R(this,T,R) result(entropy)
    class(GasSpecies),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),pointer,dimension(:)::     c
        
        if(T>this%Tc_) then
            c => this%ThermoCoefH()
        else
            c => this%ThermoCoefL()
        end if
        entropy = -0.5_rp*c(1)/T**2 - c(2)/T + c(3)*log(T) + c(9)
        entropy = entropy + sum(c(4:7)*T**[1:4]/real([1:4],kind=rp))
        entropy = R * entropy
        
    end function EntropyNasa9_R
    
    !specific entropy polynomial fitting
    elemental real(rp) function EntropyNasa9_Rc(this,T) result(entropy)
    class(GasSpecies),intent(in)::  this
    real(rp),intent(in)::           T
        entropy = this%S(T,R_c)
    end function EntropyNasa9_Rc
    
    
    
    !--refer to chemkin manual [pure species vicosity]
    elemental real(rp) function viscosity_sp(this,T) result(miu)
    class(GasSpecies),intent(in)::  this
    real(rp),intent(in)::           T
    real(rp)::                      omega, a
    
        !refer to [a fortran computer code package for the evaluation of gas-phase, multicomponent transport properties]
        a = T/this%epsilonOverk_
        omega = 1.147_rp*a**(-0.145_rp) + (a + 0.5)**(-2)
    
        !miu = (5/16)[sqrt(pi*m*k*T)/pi/sigma**2/omega]
        miu = 2.6693e-6_rp*(sqrt(this%mw_*T)/this%sigma_**2/omega)

    end function viscosity_sp
    !------------------------------------------------------------------------  
    
    
!gas chem reaction model 
!-------------------------------------------------------------------------
    pure subroutine init_sp(this, sp)
    class(GasChemReactKin),intent(out)::        this
    type(GasSpecies),dimension(:),intent(in)::  sp
        allocate(this%species_, source=sp)
    end subroutine init_sp
    
    pure subroutine init_spReact(this, sp, sm, arns, threebody, arnsb)
    class(GasChemReactKin),intent(out)::        this
    type(GasSpecies),dimension(:),intent(in)::  sp
    integer(ip),dimension(:,:,:),intent(in)::   sm
    real(rp),dimension(:,:),intent(in)::        threebody,arns
    real(rp),dimension(:,:),optional,intent(in)::arnsb
    integer(ip)::                               i,ns,nr
    
        ns = size(sp)
        nr = size(arns,2)
        
        this%ns_ = ns
        this%nr_ = nr
        
        allocate(this%species_, source=sp)
        allocate(this%stoichiometric_, source=sm)
        allocate(this%threebody_, source=threebody)
        allocate(this%arns_, source=arns)
        if(present(arnsb)) allocate(this%arnsRev_, source=arnsb)
        
    end subroutine init_spReact
    
    !--
    pure subroutine init_movealloc(this,that)
    class(GasChemReactKin),intent(out)::    this
    class(GasChemReactKin),intent(inout)::  that
        
        this%ns_ = that%ns_
        this%nr_ = that%nr_
        call move_alloc(that%species_,this%species_)
        call move_alloc(that%stoichiometric_,this%stoichiometric_)
        call move_alloc(that%threebody_,this%threebody_)
        call move_alloc(that%arns_,this%arns_)
        if(allocated(that%arnsRev_)) call move_alloc(that%arnsRev_,this%arnsRev_)
    
    end subroutine init_movealloc
    
    
    !this is a corrected version, details refer to wiki Arrhenius method
    elemental real(rp) function Arrhenius(a,b,Et,T) result(k)
    real(rp),intent(in)::                       a,b,Et,T
    
        k = a * T**b * exp(-Et/T)
        
    end function Arrhenius
    
    elemental real(rp) function Arrheniusf_i(this,i,T) result(k)
    class(GasChemReactKin),intent(in)::         this
    integer(ip),intent(in)::                    i
    real(rp),intent(in)::                       T
    
        k = this%k(this%arns_(1,i), this%arns_(2,i), this%arns_(3,i), T)
        
    end function Arrheniusf_i
    
    !production rate of concentraction[mol], Species concentraction[mol] = rho_i / MolecularWeight
    pure function ProdRate_Basic(MolCon,smf,smr,kf,kr,TbCon) result(pr)
    real(rp),dimension(:),intent(in)::          MolCon,kf,kr,TbCon
    integer(ip),dimension(:,:),intent(in)::     smf,smr
    real(rp),dimension(size(MolCon))::          pr
    integer(ip)::                               si
    
        do si=1,size(pr)
            pr(si) = sum((smr(:,si) - smf(:,si)) * ProgRate(MolCon,smf,smr,kf,kr) * TbCon)
        enddo
        
    end function ProdRate_Basic
    
    !--
    pure function ProdRate_mc(this,MolCon,T) result(pr)
    class(GasChemReactKin),intent(in)::         this
    real(rp),dimension(:),intent(in)::          MolCon  !mol/cm3
    real(rp),intent(in)::                       T
    real(rp),dimension(size(molcon))::          pr
    real(rp),dimension(:),allocatable::         kf,kr
    
        allocate(kf, source = this%kf([1:this%nr()], T))
        allocate(kr, source = this%Revk(kf, this%sp(), this%smf(), this%smr(), T))
        pr = this%molPr(Molcon, this%smf(), this%smr(), kf, kr, this%TbCon(MolCon))
        
    end function ProdRate_mc
    
    !rate of progress variables
    pure function ProgRate(MolCon,smf,smr,kf,kr)
    real(rp),dimension(:),intent(in)::          MolCon,kf,kr
    integer(ip),dimension(:,:),intent(in)::     smf,smr
    real(rp),dimension(size(kF))::              ProgRate
    
        ProgRate = prograte1(molcon,smf,kf) - prograte1(molcon,smr,kr)
        
    end function ProgRate
    
    !1 direction prograte
    pure function ProgRate1(molcon,sm,k)
    real(rp),dimension(:),intent(in)::          molcon,k
    integer(ip),dimension(:,:),intent(in)::     sm
    real(rp),dimension(size(k))::               ProgRate1
    integer(ip)::                               i
    
        do i=1,size(progRate1)  !nr
            progRate1(i) = k(i)*product(molcon**sm(i,:))
        enddo
        
    end function ProgRate1

    
    !RevRateConstant
    pure function RevRateConstant_equilibrium(rateConstant,sp,smf,smr,T) result(rrc)
    real(rp),dimension(:),intent(in)::          rateConstant
    type(GasSpecies),dimension(:),intent(in)::  sp
    integer(ip),dimension(:,:),intent(in)::     smf,smr
    real(rp),intent(in)::                       T
    real(rp),dimension(size(rateConstant))::    rrc
    
        rrc = rateConstant / EquilibriumConstant(sp,smr-smf,T)
        
    end function RevRateConstant_equilibrium

    
    !dsm = smr - smf | (nr,ns)
    !pressure should be scaled from pa -> dynes/cm2 | p_dynes = p_atm * 10
    !universe gas constant should be scaled from J/(mol K) -> ergs/(mol K) | R_ergs = R_c * 10^7
    !then Patm / R should be scalsed like p_dynes/R_ergs = P_atm / R_c * 10^(-6) 
    !refer to Chemkin theory mannual
    pure function EquilibriumConstant(sp,dsm,T) result(c)
    type(GasSpecies),dimension(:),intent(in)::  sp
    integer(ip),dimension(:,:),intent(in)::     dsm
    real(rp),intent(in)::                       T
    real(rp),dimension(size(dsm,1))::           c
    integer(ip)::                               k,ns
    real(rp),parameter::                        Teq = P_atm / R_c * 1.e-6_rp
    
        ns = size(sp)
        do k = 1,size(c) !nr
            c(k) = exp(sum(dsm(k,:) * (sp([1:ns])%S(T) - sp([1:ns])%H(T)/T)) / R_c) !Kpi
            c(k) = c(k) * (Teq/T)**(sum(dsm(k,:))) !Kci
        enddo
        
    end function EquilibriumConstant
    
    !Three body concentration = sum(molality * Three body coefficient)
    pure function ThreebodyCon(MolCon,ThreeBody) result(TbCon)
    real(rp),dimension(:),intent(in)::          MolCon ! [mol/cm^3] mole concentration
    real(rp),dimension(:,:),intent(in)::        ThreeBody
    real(rp),dimension(size(ThreeBody,1))::     TbCon
    integer(ip)::                               re,nr
        
        nr = size(ThreeBody, 1)
        do re =1,nr
            if(abs(sum(ThreeBody(re,:))) < Globaleps) then
                TbCon(re) = 1._rp !not three body reaction
            else
                TbCon(re) = sum(ThreeBody(re,:)*MolCon)
            endif
        end do
          
    end function ThreebodyCon
    
    !--
    pure function ThreebodyCon_molcon(this,MolCon) result(TbCon)
    class(GasChemReactKin),intent(in)::         this
    real(rp),dimension(:),intent(in)::          MolCon ! [mol/cm^3] mole concentration
    real(rp),dimension(:),allocatable::         TbCon
        
        allocate(TbCon,source=this%TbCon(MolCon,this%ThreeBody_))

    end function ThreebodyCon_molcon
    
    !caculate gas temperature 
    !Ri = Rc/mw [J][K-1][mol-1][mol][g-1] = [J][K-1][g-1]
    pure real(rp) function Temperature_Sp(sp,T0,e,massfrac) result(T)
    type(GasSpecies),dimension(:),intent(in)::  sp          !
    real(rp),intent(in)::                       T0,e        !e = Cv*T inner e [J/g]
    real(rp),dimension(:),intent(in)::          massfrac    !mass fraction
    integer(ip)::                               counter
    real(rp)::                                  Tt,f,fp
    real(rp),dimension(size(sp))::              cpi,hi,Ri
    
        Ri = sp%R()
        
        !newton iterative method
        ![f] is energy conservative function
        ![fp] is the derived function
        counter=0; T = T0; Tt = 0._rp
        do while(counter<=15 .and. abs(T-Tt) > GlobalEps*100._rp)
            counter = counter + 1
            Tt = T
            Cpi = sp%Cp(Tt,Ri); Hi  = sp%H(Tt,Ri)
            
            !f = rhoH - rhoRT - e = 0._rp | e is constant | J/g
            f  = sum(massfrac*Hi) - sum(massfrac*Ri)*Tt - e
            !fp = rhoCp - rhoR
            fp = sum(massfrac*Cpi) - sum(massfrac*Ri)
            !x = x0 - f(x)/f'(x)
            T = Tt - f/fp
        end do
        
    end function Temperature_Sp
    
    pure real(rp) function temperature_kin(this,T0,e,massfrac) result(T)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T0,e        !e = Cv*T inner e [J/g]
    real(rp),dimension(:),intent(in)::      massfrac    !mass fraction
    type(GasSpecies),dimension(:),pointer:: sp          
        
        sp => this%sp()
        T = this%T(sp,T0,e,massfrac)
    
    end function temperature_kin
    
    !--output j/(K*g)
    pure real(rp) function Cp_mf(this,T,massfrac) result(c)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T
    real(rp),dimension(:),intent(in)::      massfrac
    integer(ip)::                           i
        c = sum(massfrac * this%species_%Cp(T,this%species_%R()))   
    end function Cp_mf
    
    !--output j/g
    pure real(rp) function enthalpy_mf(this,T,massfrac) result(e)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T
    real(rp),dimension(:),intent(in)::      massfrac
    integer(ip)::                           i
        e = sum(massfrac * this%species_%H(T,this%species_%R()))
    end function enthalpy_mf
    
    !--output j/g
    pure real(rp) function entropy_mf(this,T,massfrac) result(e)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T
    real(rp),dimension(:),intent(in)::      massfrac
        e = sum(massfrac * this%species_%S(T,this%species_%R()))
    end function entropy_mf
    
    !--output j/(K*g)
    pure real(rp) function R_mf(this,massfrac)
    class(GasChemReactKin),intent(in)::     this
    real(rp),dimension(:),intent(in)::      massfrac
        R_mf = sum(massfrac * this%species_([1:this%ns_])%R())
    end function R_mf
    
    
    pure type(GasSpecies) function sp_i(this,i)
    class(GasChemReactKin),intent(in)::   this
    integer(ip),intent(in)::                i
        sp_i = this%species_(i)
    end function sp_i
    
    function sp_ptr(this)
    class(GasChemReactKin),target,intent(in)::this
    type(GasSpecies),dimension(:),pointer::sp_ptr 
        sp_ptr => this%species_
    end function sp_ptr
    
    function smf(this)
    class(GasChemReactKin),target,intent(in)::this
    integer(ip),dimension(:,:),pointer::    smf
        smf => this%stoichiometric_(:,:,1)
    end function smf
    
    function smr(this)
    class(GasChemReactKin),target,intent(in)::this
    integer(ip),dimension(:,:),pointer::    smr
        smr => this%stoichiometric_(:,:,2)
    end function smr
    
    pure integer(ip) function ns(this)
    class(GasChemReactKin),intent(in)::   this
        ns = this%ns_
    end function ns
    
    pure integer(ip) function nr(this)
    class(GasChemReactKin),intent(in)::   this
        nr = this%nr_
    end function nr
    
    
!---------------------------------
    pure real(rp) function miu_mix(sp,molfrac,T) result(miu)
    type(GasSpecies),dimension(:),intent(in)::  sp
    real(rp),dimension(:),intent(in)::          molfrac
    real(rp),intent(in)::                       T
    integer(ip)::                               i,j,ns
    real(rp),dimension(size(sp),size(sp))::     phi
    real(rp),dimension(size(sp))::              miuk
    
        miuk = sp%miu(T)
        ns = size(sp)
        
        do i=1,ns
            do j=1,ns
                phi(i,j) = (1._rp + sqrt(miuk(i)/miuk(j))*(sp(j)%mw_/sp(i)%mw_)**0.25_rp)**2
                phi(i,j) = phi(i,j)/sqrt(8._rp*(1._rp + sp(i)%mw_/sp(j)%mw_))
            enddo
        enddo
        
        miu = 0._rp
        do i=1,ns
            miu = miu + molfrac(i)*miuk(i)/sum(molfrac*phi(i,:))
        enddo
        
    end function miu_mix
    
end module GasChemReactKin_