#include "AUSMPWPLUSFlux.H"

scalar Foam::AUSMPWPLUSFlux::funF
(
        const scalar& p,
        const scalar& pf,
        const scalar& M
) const
{
    if(mag(M) <= 1)
    {
        return p/pf - 1;
    }
    else
    {
        return 0;
    }
}

scalar Foam::AUSMPWPLUSFlux::funw
(
    const scalar& p1,
    const scalar& p2
) const
{
    return 1 - pow(min(p1/p2, p2/p1), 3);
}

void Foam::AUSMPWPLUSFlux::evaluateFlux
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
    const scalar cSf_pos = c_pos*magSf;
    const scalar cSf_neg = c_neg*magSf;

    aByU = (U_pos + U_neg)/2.0;
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);

    const scalar h_pos = e_pos + 0.5 * magSqr(U_pos) + p_pos/rho_pos;
    const scalar h_neg = e_neg + 0.5 * magSqr(U_neg) + p_neg/rho_neg;

    const vector unit_normal = Sf / magSf;
    const scalar u_l = U_pos & unit_normal;
    const scalar u_r = U_neg & unit_normal;

    const scalar v_l2 = magSqr(U_pos) - sqr(u_l);
    const scalar v_r2 = magSqr(U_neg) - sqr(u_r);

    const scalar h_n = 0.5*(h_pos - 0.5*v_l2 + h_neg - 0.5*v_r2);

    const scalar cstar_pos = Foam::sqrt(max(((2*(gamma_pos - 1)/(gamma_pos + 1))*h_n), SMALL));
    const scalar cstar_neg = Foam::sqrt(max(((2*(gamma_neg - 1)/(gamma_neg + 1))*h_n), SMALL));

    const scalar c_f = (0.5*(u_l + u_r)) > 0 ? sqr(cstar_pos)/max(u_l, cstar_pos) : sqr(cstar_neg)/max(u_r, cstar_neg);
    //const scalar gamma = 0.5*(((p_pos/(rho_pos*e_pos)) + 1) + ((p_neg/(rho_neg*e_neg)) + 1));
    //const scalar gamma = 0.5*(gamma_pos + gamma_neg);
    //const scalar c_star = sqrt(max(2*(gamma - 1)/(gamma + 1)*h_n, SMALL));

    //const scalar c_f = (0.5*(u_l + u_r)) > 0 ? sqr(c_star)/max(u_l, c_star) : sqr(c_star)/max(u_r, c_star);

    //const scalar cbar_pos = sqr(cstar_pos)/max(u_l, cstar_pos);
    //const scalar cbar_neg = sqr(cstar_neg)/max(u_r, cstar_neg);

    //const scalar c_f = min(cbar_pos, cbar_neg);

    //const scalar c_f = 0.5*(c_pos + c_neg);
    //const scalar c_f = Foam::sqrt(c_pos * c_neg);

    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + c_f * magSf;
    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);
    //amaxSf = max(mag(phiv_pos), mag(phiv_neg));

    const scalar Ma_l = u_l/c_f;
    const scalar Ma_r = u_r/c_f;

    const scalar Ma_lmag = mag(Ma_l);
    const scalar Ma_rmag = mag(Ma_r);

    const scalar Ma_pl = ((Ma_lmag >= 1) ? (0.5*(Ma_l + Ma_lmag)) : (0.25*sqr(Ma_l + 1)));
    const scalar Ma_nr = ((Ma_rmag >= 1) ? (0.5*(Ma_r - Ma_rmag)) : (-0.25*sqr(Ma_r - 1)));

    const scalar p_pl = ((Ma_lmag >= 1) ? (0.5*(Ma_l + Ma_lmag)/Ma_l) : (0.25*sqr(Ma_l + 1)*(2 - Ma_l)));
    const scalar p_nr = ((Ma_rmag >= 1) ? (0.5*(Ma_r - Ma_rmag)/Ma_r) : (0.25*sqr(Ma_r - 1)*(2 + Ma_r)));

    const scalar Ma_f = Ma_pl + Ma_nr;
    const scalar p_f = p_pl*p_pos + p_nr*p_neg;

    //const scalar u_f = c_f * Ma_f;
    //const scalar u_fmag = c_f * mag(Ma_f);

    const scalar Fl = funF(p_pos, p_f, Ma_l);
    const scalar Fr = funF(p_neg, p_f, Ma_r);

    const scalar w = funw(p_pos, p_neg);

    const scalar Mbar_pl = Ma_f >= 0 ? Ma_pl + Ma_nr*((1 - w)*(1 + Fr) - Fl) : w*(1 + Fl)*Ma_pl;
    const scalar Mbar_nr = Ma_f >= 0 ? w*(1 + Fr)*Ma_nr : Ma_nr + Ma_pl*((1 - w)*(1 + Fl) - Fr);

    phi = (Mbar_pl*c_f*rho_pos + Mbar_nr*c_f*rho_neg)*magSf;
    phiUp = (Mbar_pl*c_f*rhoU_pos + Mbar_nr*c_f*rhoU_neg)*magSf + p_f * Sf;
    phiEp = (Mbar_pl*c_f*rho_pos*h_pos + Mbar_nr*c_f*rho_neg*h_neg)*magSf;
}

