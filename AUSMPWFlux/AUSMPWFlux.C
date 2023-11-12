#include "AUSMPWFlux.H"

scalar Foam::AUSMPWFlux::funPl
(
        const scalar& p1,
        const scalar& p2
) const
{
    const scalar idt = min(p1/p2, p2/p1);
    if(idt >= 3/4 && idt < 1)
    {
        return 4*min(p1/p2, p2/p1) - 3;
    }
    else if(idt >= 0 && idt < 3/4)
    {
        return 0;
    }

}

scalar Foam::AUSMPWFlux::funFl
(
        const scalar& pl,
        const scalar& pr,
        const scalar& pf,
        const scalar& M,
        const scalar& cf,
        const scalar& magU
) const
{
    if(mag(M) <= 1)
    {
        return (pl/pf - 1)*funPl(pl, pr)*0.25*sqr(M + 1)*min(1, pow(magU/cf, 0.25));
    }
    else
    {
        return 0;
    }
}

scalar Foam::AUSMPWFlux::funFr
(
        const scalar& pl,
        const scalar& pr,
        const scalar& pf,
        const scalar& M,
        const scalar& cf,
        const scalar& magU
) const
{
    if(mag(M) <= 1)
    {
        return (pr/pf - 1)*funPl(pr, pl)*0.25*sqr(M - 1)*min(1, pow(magU/cf, 0.25));
    }
    else
    {
        return 0;
    }
}

scalar Foam::AUSMPWFlux::funw
(
    const scalar& p1,
    const scalar& p2
) const
{
    return 1 - pow(min(p1/p2, p2/p1), 3);
}

scalar Foam::AUSMPWFlux::funPws
(
    const scalar& p1,
    const scalar& p2,
    const scalar& x,
    const scalar& y
) const
{
    return (1 - funw(p1, p2))*x + funw(p1, p2)*y;
}

vector Foam::AUSMPWFlux::funPwv
(
    const scalar& p1,
    const scalar& p2,
    const vector& x,
    const vector& y
) const
{
    return (1 - funw(p1, p2))*x + funw(p1, p2)*y;
}

void Foam::AUSMPWFlux::evaluateFlux
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

    //const scalar c_pos = Foam::sqrt(max(SMALL, gamma_pos * p_pos / max(rho_pos, eps)));
    //const scalar c_neg = Foam::sqrt(max(SMALL, gamma_neg * p_neg / max(rho_neg, eps)));

    const scalar cSf_pos = c_pos*magSf;
    const scalar cSf_neg = c_neg*magSf;

    aByU = (U_pos + U_neg)/2.0;
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);

    const scalar h_pos = e_pos + 0.5 * magSqr(U_pos) + p_pos/rho_pos;
    const scalar h_neg = e_neg + 0.5 * magSqr(U_neg) + p_neg/rho_neg;

    const scalar cstar_pos = sqrt(2*(gamma_pos - 1)/(gamma_pos + 1)*h_pos);
    const scalar cstar_neg = sqrt(2*(gamma_neg - 1)/(gamma_neg + 1)*h_neg);

    const vector unit_normal = Sf / magSf;
    const scalar u_l = U_pos & unit_normal;
    const scalar u_r = U_neg & unit_normal;
/*
    const scalar v_l = sqrt(magSqr(U_pos) - sqr(u_l));
    const scalar v_r = sqrt(magSqr(U_neg) - sqr(u_r));

    const scalar h_n = 0.5*(h_pos - 0.5*sqr(v_l) + h_neg - 0.5*sqr(v_r));

    const scalar gamma = 0.5*(gamma_pos + gamma_neg);
    const scalar c_star = sqrt(2*(gamma - 1)/(gamma + 1)*h_n);

    const scalar c_f = (0.5*(u_l + u_r)) > 0 ? sqr(c_star)/max(u_l, c_star) : sqr(c_star)/max(u_r, c_star);
*/
    const scalar cbar_pos = sqr(cstar_pos)/max(u_l, cstar_pos);
    const scalar cbar_neg = sqr(cstar_neg)/max(u_r, cstar_neg);

    const scalar c_f = min(cbar_pos, cbar_neg);

    //const scalar c_f = 0.5*(c_pos + c_neg);
    //const scalar c_f = Foam::sqrt(c_pos * c_neg);

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

    //const scalar u_f = c_f * Ma_f;
    //const scalar u_fmag = c_f * mag(Ma_f);

    const scalar Fl = funFl(p_pos, p_neg, p_f, Ma_l, c_f, mag(U_pos));
    const scalar Fr = funFr(p_pos, p_neg, p_f, Ma_r, c_f, mag(U_neg));

     scalar pws = 0;
     vector pwv = Foam::vector::zero;

    if(Ma_f >= 0)
    {
        pws = funPws(p_pos, p_neg, rho_pos, rho_neg);
        pwv = funPwv(p_pos, p_neg, rhoU_pos, rhoU_neg);

        phi = ((1 + Fl)*Ma_pl*c_f*rho_pos + (1 + Fr)*Ma_nr*c_f*pws)*magSf;
        phiUp = ((1 + Fl)*Ma_pl*c_f*rhoU_pos + (1 + Fr)*Ma_nr*c_f*pwv)*magSf + p_f * Sf;
        phiEp = ((1 + Fl)*Ma_pl*c_f*rho_pos*h_pos + (1 + Fr)*Ma_nr*c_f*pws*h_pos)*magSf;
    }
    else
    {
        pws = funPws(p_pos, p_neg, rho_neg, rho_pos);
        pwv = funPwv(p_pos, p_neg, rhoU_neg, rhoU_pos);

        phi = ((1 + Fl)*Ma_pl*c_f*pws + (1 + Fr)*Ma_nr*c_f*rho_neg)*magSf;
        phiUp = ((1 + Fl)*Ma_pl*c_f*pwv + (1 + Fr)*Ma_nr*c_f*rhoU_neg)*magSf + p_f * Sf;
        phiEp = ((1 + Fl)*Ma_pl*c_f*pws*h_neg + (1 + Fr)*Ma_nr*c_f*rho_neg*h_neg)*magSf;
    }
}

