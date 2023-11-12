#include "fvCFD.H"

namespace Foam
{



class AUSMPLUSFlux
{
public:

    // Constructors

        //- Constructor null
        AUSMPLUSFlux()
        {}


    // Destructor - default


    // Member Functions

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

