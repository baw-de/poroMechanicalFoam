/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    SANISAND

Description
    Elastic-plastic mechanical law where the elasticity is given by Hooke's law
    and the plasticity by Mohr Coulomb.

    Note: to include poro-pressure (poro-elasticity) effects, use the
    poroLinearElasticMohrCoulomb law.

    More reading at:

    T.Tang, O.Hededal & P.Cardiff (2014). On finite volume method implementation
    of poro-elasto-plasticity soil model. International Journal for
    Numerical and Analytical Methods in Geomechanics, 10.1002/nag.2361.

    T.Tang & O.Hededal (2014). Simulation of pore pressure accumulation under
    cyclic loading using finite volume method. Proceedings of the 8th European
    Conference on Numerical Methods in Geotechnical Engineering (NUMGE14),
    Volume 2, Pages 1301-1306.

SourceFiles
    SANISAND.C

Author
    Tian Tang, DTU. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.
    Denis Maier, BAW: brought into solids4Foam format.

\*---------------------------------------------------------------------------*/

#ifndef SANISAND_H
#define SANISAND_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SANISAND Declaration
\*---------------------------------------------------------------------------*/

class SANISAND
:
    public mechanicalLaw
{
    // Private data

    // Properties

        //- Density
    dimensionedScalar rho_;
    dimensionedScalar e0_;
    dimensionedScalar K0_;
    dimensionedScalar G0_;
    dimensionedScalar lambda0_;
    dimensionedScalar e_cref_;
    dimensionedScalar Pref_;
    dimensionedScalar a_;
    dimensionedScalar Mc_;
    dimensionedScalar Me_;
    dimensionedScalar LAMBDA_;
    dimensionedScalar A0_;
    dimensionedScalar Cf_;
    dimensionedScalar Fmax_;
    dimensionedScalar h0_;
    dimensionedScalar m0_;
    dimensionedScalar cm_;
    dimensionedScalar k_cb_;
    dimensionedScalar k_eb_;
    dimensionedScalar k_cd_;
    dimensionedScalar k_ed_;
    dimensionedScalar c_;
    dimensionedScalar c_b_;
    dimensionedScalar c_d_;
    scalar tol_f_;
    scalar plast_relax_;
    scalar nCorr_;

    // Fields

    volSymmTensorField DSigma_;
    volScalarField P_;
    volSymmTensorField s_;
    volSymmTensorField r_;
    volSymmTensorField alpha_;
    volSymmTensorField F_;
    volScalarField m_;
    volSymmTensorField DEpsilon_;
    volScalarField DEpsilonV_;
    volSymmTensorField DEpsilonD_;
    volScalarField epsilonV_;
    volSymmTensorField epsilonD_;

    //- Incremental change of plastic strain
    volSymmTensorField DEpsilonP_;

    volScalarField DEpsilonPV_;
    volSymmTensorField DEpsilonPD_;

    //- Total accumulated plastic strain
    volSymmTensorField epsilonP_;

    //- void ratio
    volScalarField e_;
    volScalarField psi_;
    //- Bulk modulus
    volScalarField K_;
    //- Second Lame parameter (shear modulus)
    volScalarField mu_;

    //- Active yielding flag
    //     1.0 for active yielding
    //     0.0 otherwise
    volScalarField activeYield_;

    // Private Member Functions

        //- Calculate the stress
    void calculateStress(
        symmTensor &stress,
        const scalar &dEpsV,
        const symmTensor &dEpsD,
        const scalar &poreRatio,
        const scalar &POld,
        const symmTensor &sOld,
        const symmTensor &alphaOld,
        const symmTensor &FOld,
        const scalar &mOld,
        scalar &dEpsPV,
        symmTensor &dEpsPD,
        scalar &PNew,
        symmTensor &sNew,
        symmTensor &alphaParam,
        scalar &mParam,
        symmTensor &FParam,
        scalar &KParam,
        scalar &GParam,
        scalar &psiParam,
        scalar &yieldFlag) const;

    //- Disallow default bitwise copy construct
    SANISAND(const SANISAND &);

    //- Disallow default bitwise assignment
    void operator=(const SANISAND &);

public:

    //- Runtime type information
    TypeName("SANISAND");

    // Static data members


    // Constructors

        //- Construct from dictionary
        SANISAND
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~SANISAND();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar& mu() const;

        //- Provide access to lambda for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar& lambda() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
