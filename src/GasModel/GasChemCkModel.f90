module GasChemReactKinCk_
use constants, only: ip, rp, lp, cl, R_c, P_atm, GlobalEps
use stringOpsLib
use chemkinWrapperLib
use GasSpecies_
use GasChemReactKin_
implicit none

    private
    public:: GasChemReactKinCk

    !----
    type,extends(GasChemReactKin):: GasChemReactKinCk

		type(ckMech)::			ck_

    contains

		generic::				init => init_ck
		procedure::				init_ck

    end type GasChemReactKinCk

!--
contains

    subroutine init_ck(this, fn)
	class(GasChemReactKinCk),intent(out)::  this
	character(cl),intent(in)::				fn
	
		call this%ck_%ckinit(trim(fn))
	
	end subroutine init_ck

end module GasChemReactKinCk_
