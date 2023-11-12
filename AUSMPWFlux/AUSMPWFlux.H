#include "fvCFD.H"

namespace Foam
{



class AUSMPWFlux
{
public:

    // Constructors

        //- Constructor null
        AUSMPWFlux()
        {}


    // Destructor - default


    // Member Functions
        scalar funPl
        (
            const scalar& p1,
            const scalar& p2
        ) const;

        scalar funFl
        (
                const scalar& pl,
                const scalar& pr,
                const scalar& pf,
                const scalar& Ma,
                const scalar& cf,
                const scalar& magU
        ) const;

        scalar funFr
        (
                const scalar& pl,
                const scalar& pr,
                const scalar& pf,
                const scalar& Ma,
                const scalar& cf,
                const scalar& magU
        ) const;

        scalar funw
        (
            const scalar& p1,
            const scalar& p2
        ) const;

        scalar funPws
        (
            const scalar& p1,
            const scalar& p2,
            const scalar& x,
            const scalar& y
        ) const;

        vector funPwv
        (
            const scalar& p1,
            const scalar& p2,
            const vector& x,
            const vector& y
        ) const;

        //- Update fluxes based on current left and right states
       void evaluateFlux
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
	) const;
};

}

