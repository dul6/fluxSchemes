#include "KNPFlux.H"

void Foam::KNPFlux::evaluateFlux
(
        const vector& Sf,
        const scalar& magSf,
        const scalar& c_pos,
        const scalar& c_neg,
        const scalar& rho_pos,
        const scalar& rho_neg,
        const vector& U_pos,
        const vector& U_neg,
        const vector& rhoU_pos,
        const vector& rhoU_neg,
        const scalar& p_pos,
        const scalar& p_neg,
        const scalar& e_pos,
        const scalar& e_neg,
        const scalar& gamma_pos,
        const scalar& gamma_neg,
        vector& aByU,
        scalar& amaxSf,
        scalar& phi,
        vector& phiUp,
        scalar& phiEp
) const
{
    scalar phiv_pos = U_pos & Sf;
    scalar phiv_neg = U_neg & Sf;

    const scalar cSf_pos = c_pos * magSf;
    const scalar cSf_neg = c_neg * magSf;

    const scalar ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), 0.0);
    
    const scalar am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), 0.0);
    
    const scalar a_pos = ap/(ap - am);
   
    const scalar aSf = am*a_pos;

    scalar a_neg = 1.0 - a_pos;
    
    aByU = (a_pos*U_pos + a_neg*U_neg);
    phiv_pos *= a_pos;
    phiv_neg *= a_neg; 

    const scalar aphiv_pos = phiv_pos - aSf; 
    const scalar aphiv_neg = phiv_neg + aSf; 

    // Reuse amaxSf for the maximum positive and negative fluxes
    // estimated by the central scheme
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));
    
    // --- Fluxes
    phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;
    phiUp =
        (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
      + (a_pos*p_pos + a_neg*p_neg)*Sf;

    phiEp =
        aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
      + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
      + aSf*p_pos - aSf*p_neg;

}
