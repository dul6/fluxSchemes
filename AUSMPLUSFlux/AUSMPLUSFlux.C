#include "AUSMPLUSFlux.H"

void Foam::AUSMPLUSFlux::evaluateFlux
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
    const scalar alpha = 3/16;
    const scalar beta = 1/8;

    const scalar cSf_pos = c_pos*magSf;
    const scalar cSf_neg = c_neg*magSf;

    aByU = (U_pos + U_neg)/2.0;
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);

    const scalar rhoH_pos = rho_pos*(e_pos + 0.5 * magSqr(U_pos)) + p_pos;
    const scalar rhoH_neg = rho_neg*(e_neg + 0.5 * magSqr(U_neg)) + p_neg;

    const vector unit_normal = Sf / magSf;
    //const scalar u_l = U_pos & unit_normal;
    //const scalar u_r = U_neg & unit_normal;

    const scalar u_l = U_pos.x();
    const scalar u_r = U_neg.x();

    //const scalar c_f = 0.5*(c_pos + c_neg);
    const scalar c_f = Foam::sqrt(c_pos * c_neg);

    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + c_f * magSf;
    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);
    //amaxSf = max(mag(phiv_pos), mag(phiv_neg));

    const scalar Ma_l = u_l/c_f;
    const scalar Ma_r = u_r/c_f;

    const scalar Ma_lmag = mag(Ma_l);
    const scalar Ma_rmag = mag(Ma_r);

    const scalar Ma_pl = ((Ma_lmag >= 1) ? (0.5*(Ma_l + Ma_lmag)) : (0.25*sqr(Ma_l + 1) + beta*sqr(sqr(Ma_l) - 1)));
    const scalar Ma_nr = ((Ma_rmag >= 1) ? (0.5*(Ma_r - Ma_rmag)) : (-0.25*sqr(Ma_r - 1) - beta*sqr(sqr(Ma_r) - 1)));

    const scalar p_pl = ((Ma_lmag >= 1) ? (0.5*(Ma_l + Ma_lmag)/Ma_l) : (0.25*sqr(Ma_l + 1)*(2 - Ma_l) + alpha*Ma_l*sqr(sqr(Ma_l) - 1)));
    const scalar p_nr = ((Ma_rmag >= 1) ? (0.5*(Ma_r - Ma_rmag)/Ma_r) : (0.25*sqr(Ma_r - 1)*(2 + Ma_r) - alpha*Ma_r*sqr(sqr(Ma_r) - 1)));

    const scalar Ma_f = Ma_pl + Ma_nr;
    const scalar p_f = p_pl*p_pos + p_nr*p_neg;

    const scalar u_f = c_f * Ma_f;
    const scalar u_fmag = c_f * mag(Ma_f);

    phi = 0.5*(u_f*(rho_pos + rho_neg) - u_fmag*(rho_neg - rho_pos))*magSf;
    phiUp = 0.5*(u_f*(rhoU_pos + rhoU_neg) - u_fmag*(rhoU_neg - rhoU_pos))*magSf + p_f*Sf;
    phiEp = 0.5*(u_f*(rhoH_pos + rhoH_neg) - u_fmag*(rhoH_neg - rhoH_pos))*magSf;

}

