module GasChemReactKin_
use constants, only: ip, rp, lp, cl, R_c, P_atm, GlobalEps
use stringOpsLib
use chemkinWrapperLib
use GasSpecies_
implicit none

    private
    public:: GasChemReactKin
    !--
    public:: mu_mix
	public:: molcon

    !----
    type GasChemReactKin

        private
        integer(ip)::               ns_ = 0
		integer(ip)::				nr_ = 0

        type(GasSpecies),dimension(:),&
        allocatable::               species_

		!(nr,ns,ix) [ix=1 forward][ix=2 backward]
        integer(ip),dimension(:,:,:),&
        allocatable::               stoichiometric_

		!(nr,ns)coef for ThreebodyCon
        real(rp),dimension(:,:),&
		allocatable::   			threebody_
		!(3,nr) coef for arrhenius
        real(rp),dimension(:,:),&
		allocatable::   			arns_

		!(3,nr) coef for Rev rateconstant
        real(rp),dimension(:,:),&
		allocatable::   			arnsRev_

    contains

        generic::                   init => init_sp, init_spReact
        procedure::                 init_sp
        procedure::                 init_spReact

		!--
        procedure::                 init_movealloc

        !Chemistry production of basic model, see ChemKin of Ansys
        !mol generation rate[mol/cm3/s]
        generic::                   molPr => ProdRate_Basic, ProdRate_mc
        procedure,nopass::          ProdRate_Basic
        procedure::                 ProdRate_mc

        !mol/cm3/s
        procedure,nopass::          ProgRate
        procedure,nopass::          ProgRate1

        !calculate rate constant => k = a T^b \exp{-\frac{Ea}{T}}
        generic::                   k => Arrhenius
        procedure,nopass::          Arrhenius

		!--
		procedure,nopass::          unitConvertA

		!--
        procedure::                 kf => Arrheniusf_i
        procedure::                 Arrheniusf_i

        !calculate Rev Rate Constant based on Equilibrium assumption
        generic::                   Revk => RevRateConstant_equilibrium
        procedure,nopass::          RevRateConstant_equilibrium
        !procedure,nopass::          RevRateConstant_Saha

        !calculate Equilibrium Constant based on experience formula
        procedure,nopass::          EquilibriumConstant

        !--
        generic::                   TbCon => ThreebodyCon, ThreebodyCon_molcon
        procedure,nopass::          ThreebodyCon
        procedure::                 ThreebodyCon_molcon

        !--
        generic::                   T => temperature_sp, temperature_kin
        procedure,nopass::          temperature_sp
        procedure::                 temperature_kin

        !--
		generic::					Cp => Cp_mf
        procedure::                 Cp_mf     		!j/(K*g)
		!procedure::					cp_mf_ck

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

!--
contains



    !------------------------------------------------------------------------
	!gas chem reaction model
	!-------------------------------------------------------------------------
    pure subroutine init_sp(this, sp)
    class(GasChemReactKin),intent(out)::			this
    type(GasSpecies),dimension(:),intent(in)::		sp

        allocate(this%species_, source = sp)

    end subroutine init_sp

	!--
    pure subroutine init_spReact(this, sp, sm, arns, threebody, arnsb)
    class(GasChemReactKin),intent(out)::			this
    type(GasSpecies),dimension(:),intent(in)::		sp
    integer(ip),dimension(:,:,:),intent(in)::		sm
    real(rp),dimension(:,:),intent(in)::			threebody,arns
    real(rp),dimension(:,:),optional,intent(in)::	arnsb
    integer(ip)::									i,ns,nr

        ns = size(sp)
        nr = size(arns, 2)

        this%ns_ = ns
        this%nr_ = nr

        allocate(this%species_, source = sp)
        allocate(this%stoichiometric_, source = sm)
        allocate(this%threebody_, source = threebody)
        allocate(this%arns_, source = arns)
        if(present(arnsb)) allocate(this%arnsRev_, source = arnsb)

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
        if(allocated(that%arnsRev_)) call move_alloc(that%arnsRev_, this%arnsRev_)

    end subroutine init_movealloc


    !this is a corrected version, details refer to wiki Arrhenius method
    elemental real(rp) function Arrhenius(a,b,Et,T) result(k)
    real(rp),intent(in)::	a,b,Et,T

        k = a*T**b*exp(-Et/T)

    end function Arrhenius

	!cgs=(mol, s, K, cm) | si=(mol, s, K, m) | if kmolidx=true, si = (kmol, s, K, m)
	!d[C_i]/dt = a*T**b*exp(-Et/T) \prod [C_j]^eidx_j
	!sum_expidx = sum(eidx_j)
	!unit of A is (mol/cm^3)^{1 - sum_expidx}/s/K^b -> (mol/m^3)^{1 - sum_expidx}/s/K^b
	!                                           or (kmol/m^3)^{1 - sum_expidx}/s/K^b
	!if(thirdidx) | {1 - sum_expidx} -> {- sum_expidx}
	!if(cvdir) si->cgs else cgs->si endif
	real(rp) function unitConvertA(a, sum_expidx, cvdir, kmolidx, thirdidx) result(ac)
	real(rp),intent(in)::       a, sum_expidx
	logical(lp),intent(in)::    cvdir, kmolidx, thirdidx
	real(rp)::		            c

	    !--
	    if(kmolidx) then
		    c = 3._rp
		else
		    c = 6._rp
		endif

		!--
		if(cvdir) c = -c

		!--
	    if(thirdidx) then
			ac = 10._rp**(- c*sum_expidx)*a
		else
			ac = 10._rp**(c - c*sum_expidx)*a
		endif

	end function unitConvertA

	!--
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
            pr(si) = sum((smr(:,si) - smf(:,si))*ProgRate(MolCon, smf, smr, kf, kr)*TbCon)
        enddo

    end function ProdRate_Basic

    !--
    function ProdRate_mc(this,MolCon,T) result(pr)
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
    function RevRateConstant_equilibrium(rateConstant,sp,smf,smr,T) result(rrc)
    real(rp),dimension(:),intent(in)::          rateConstant
    type(GasSpecies),dimension(:),intent(in)::  sp
    integer(ip),dimension(:,:),intent(in)::     smf,smr
    real(rp),intent(in)::                       T
    real(rp),dimension(size(rateConstant))::    rrc

        rrc = rateConstant/EquilibriumConstant(sp, smr-smf, T)

    end function RevRateConstant_equilibrium


    !dsm = smr - smf | (nr,ns)
    !pressure should be scaled from pa -> dynes/cm2 | p_dynes = p_atm * 10
    !universe gas constant should be scaled from J/(mol K) -> ergs/(mol K) | R_ergs = R_c*10^7
    !then Patm/R should be scalsed like p_dynes/R_ergs = P_atm/R_c*10^(-6)
    !refer to Chemkin theory mannual
    function EquilibriumConstant(sp,dsm,T) result(c)
    type(GasSpecies),dimension(:),intent(in)::  sp
    integer(ip),dimension(:,:),intent(in)::     dsm
    real(rp),intent(in)::                       T
    real(rp),dimension(size(dsm,1))::           c
    integer(ip)::                               k,ns
    real(rp),parameter::                        Teq = P_atm/R_c*1.e-6_rp

        ns = size(sp)
        do k = 1,size(c) !nr
            c(k) = exp(sum(dsm(k,:)*(sp([1:ns])%S(T) - sp([1:ns])%H(T)/T))/R_c) !Kpi
            c(k) = c(k)*(Teq/T)**(sum(dsm(k,:))) !Kci
        enddo

    end function EquilibriumConstant

    !Three body concentration = sum(molality * Three body coefficient)
    function ThreebodyCon(MolCon,ThreeBody) result(TbCon)
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
    function ThreebodyCon_molcon(this,MolCon) result(TbCon)
    class(GasChemReactKin),intent(in)::         this
    real(rp),dimension(:),intent(in)::          MolCon ! [mol/cm^3] mole concentration
    real(rp),dimension(:),allocatable::         TbCon

        allocate(TbCon,source=this%TbCon(MolCon,this%ThreeBody_))

    end function ThreebodyCon_molcon

    !caculate gas temperature
    !Ri = Rc/mw [J][K-1][mol-1][mol][g-1] = [J][K-1][g-1]
    real(rp) function Temperature_Sp(sp, T0, e, massfrac) result(T)
    type(GasSpecies),dimension(:),intent(in)::  sp          !
    real(rp),intent(in)::                       T0, e       !e = Cv*T inner e [J/g]
    real(rp),dimension(:),intent(in)::          massfrac    !mass fraction
    integer(ip)::                               counter
    real(rp)::                                  Tt,f,fp
    real(rp),dimension(size(sp))::              cpi,hi,Ri

        Ri = sp%R()
        !newton iterative method
        ![f] is energy conservative function
        ![fp] is the derived function
        counter=0
		T = T0
		Tt = 0._rp
        do while(counter<=15 .and. abs(T-Tt) > GlobalEps*100._rp)

            counter = counter + 1
            Tt = T
            Cpi = sp%Cp(Tt, Ri)
			Hi = sp%H(Tt, Ri)

            !f = rhoH - rhoRT - e = 0._rp | e is constant | J/g
            f = sum(massfrac*Hi) - sum(massfrac*Ri)*Tt - e
            !fp = rhoCp - rhoR
            fp = sum(massfrac*Cpi) - sum(massfrac*Ri)
            !x = x0 - f(x)/f'(x)
            T = Tt - f/fp

        end do

    end function Temperature_Sp

	!--
    real(rp) function temperature_kin(this, T0, e, massfrac) result(T)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T0, e       !e = Cv*T inner e [J/g]
    real(rp),dimension(:),intent(in)::      massfrac    !mass fraction
    type(GasSpecies),dimension(:),pointer:: sp

        sp => this%sp()
        T = this%T(sp, T0, e, massfrac)

    end function temperature_kin

    !--output j/(K*g)
    real(rp) function Cp_mf(this,T,massfrac) result(c)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T
    real(rp),dimension(:),intent(in)::      massfrac
    integer(ip)::                           i

        c = sum(massfrac*this%species_%Cp(T, this%species_%R()))

    end function Cp_mf


    !--output j/g
    real(rp) function enthalpy_mf(this,T,massfrac) result(e)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T
    real(rp),dimension(:),intent(in)::      massfrac
    integer(ip)::                           i

        e = sum(massfrac*this%species_%H(T, this%species_%R()))

    end function enthalpy_mf

    !--output j/g
    real(rp) function entropy_mf(this,T,massfrac) result(e)
    class(GasChemReactKin),intent(in)::     this
    real(rp),intent(in)::                   T
    real(rp),dimension(:),intent(in)::      massfrac

        e = sum(massfrac*this%species_%S(T, this%species_%R()))

    end function entropy_mf

    !--output j/(K*g)
	!R = Rc/(\sum X_i M_i) = \sum(Y_i/M_i) Rc = \sum Y_i R_i
    pure real(rp) function R_mf(this, massfrac)
    class(GasChemReactKin),intent(in)::     this
    real(rp),dimension(:),intent(in)::      massfrac

        R_mf = sum(massfrac*this%species_([1:this%ns_])%R())

    end function R_mf


    pure type(GasSpecies) function sp_i(this,i)
    class(GasChemReactKin),intent(in)::   	this
    integer(ip),intent(in)::                i

        sp_i = this%species_(i)

    end function sp_i

    function sp_ptr(this)
    class(GasChemReactKin),target,intent(in)::	this
    type(GasSpecies),dimension(:),pointer::		sp_ptr

        sp_ptr => this%species_

    end function sp_ptr

    function smf(this)
    class(GasChemReactKin),target,intent(in)::	this
    integer(ip),dimension(:,:),pointer::    	smf

        smf => this%stoichiometric_(:,:,1)

    end function smf

    function smr(this)
    class(GasChemReactKin),target,intent(in)::	this
    integer(ip),dimension(:,:),pointer::    	smr

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

    !--
    pure real(rp) function mu_mix(sp, molfrac, T) result(mu)
    type(GasSpecies),dimension(:),intent(in)::  sp
    real(rp),dimension(:),intent(in)::          molfrac
    real(rp),intent(in)::                       T
    integer(ip)::                               i,j,ns
    real(rp),dimension(size(sp),size(sp))::     phi
    real(rp),dimension(size(sp))::              muk

        muk = sp%mu(T)
        ns = size(sp)

        do i=1,ns
            do j=1,ns
                phi(i,j) = (1._rp + sqrt(muk(i)/muk(j))*(sp(j)%mw()/sp(i)%mw())**0.25_rp)**2
                phi(i,j) = phi(i,j)/sqrt(8._rp*(1._rp + sp(i)%mw()/sp(j)%mw()))
            enddo
        enddo

        mu = 0._rp
        do i=1,ns
            mu = mu + molfrac(i)*muk(i)/sum(molfrac*phi(i,:))
        enddo

    end function mu_mix

	!output[mol/cm3]
	pure subroutine molcon(sp, massfrac, rho, molc)
	type(GasSpecies),dimension(:),intent(in)::  sp
	real(rp),dimension(:),intent(in)::			massfrac
	real(rp),intent(in)::						rho
	real(rp),dimension(:),intent(out)::			molc

		![kg/m3] -> [g/cm3]
		molc = rho*massfrac*1.e-3_rp
		![g/cm3]/[g/mol] = [mol/cm3]
		molc = molc/sp%mw()

	end subroutine molcon

end module GasChemReactKin_