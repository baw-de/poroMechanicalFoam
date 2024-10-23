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
    linearElasticMohrCoulombPlasticDilationCutoff

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
    linearElasticMohrCoulombPlasticDilationCutoff.C

Author
    Tian Tang, DTU. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.
    Denis Maier, BAW. Added dilation cut-off.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMohrCoulombPlasticDilationCutoff_H
#define linearElasticMohrCoulombPlasticDilationCutoff_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                         Class linearElasticMohrCoulombPlasticDilationCutoff Declaration
\*---------------------------------------------------------------------------*/

    class linearElasticMohrCoulombPlasticDilationCutoff
        : public mechanicalLaw
    {
        // Private data

        //- Density
        const dimensionedScalar rho_;

        //- Modulus of elasticity
        const dimensionedScalar E_;

        //- Poisson's ratio
        const dimensionedScalar nu_;

        //- First Lame parameter
        const dimensionedScalar lambda_;

        //- Second Lame parameter (shear modulus)
        const dimensionedScalar mu_;

        //- Bulk modulus
        const dimensionedScalar K_;

        //- Friction angle
        const dimensionedScalar varPhi_;

        //- Cohesion
        const dimensionedScalar c_;

        //- Dilation angle
        const dimensionedScalar varPsi0_;

        //- initial void ratio
        const dimensionedScalar e0_;

        //- critical void ratio
        const dimensionedScalar eCrit_;

        //- current void ratio
        volScalarField e_;

        //- Derived plasticity parameter
        const scalar k_;

        //- Derived plasticity parameter
        const scalar m0_;

        //- Derived plasticity parameter
        const vector a_;

        //- Derived plasticity parameter
        const vector b0_;

        //- Derived plasticity parameter
        const symmTensor C_;

        //- Store inverse of C for efficiency
        const symmTensor invC_;

        //- Derived plasticity parameter
        const vector rp0_;

        //- Derived plasticity parameter
        const vector r_lf1_;

        //- Derived plasticity parameter
        const vector r_lf2_;

        //- Derived plasticity parameter
        const vector r_lg10_;

        //- Derived plasticity parameter
        const vector r_lg20_;

        //- Derived plasticity parameter
        const vector sigma_a_;

        //- Effective stress
        volSymmTensorField Dsigma_;

        //- Effective stress surface field
        surfaceSymmTensorField Dsigmaf_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Incremental change of plastic strain surface field
        surfaceSymmTensorField DEpsilonPf_;

        //- Total accumulated plastic strain
        volSymmTensorField epsilonP_;

        //- Equivalent plastic strain
        volScalarField epsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;
        volScalarField f_;
        volScalarField EigenError_;

        // Private Member Functions

        //- Calculate the stress
        void calculateStress(symmTensor &sigma, scalar &activeYield, scalar &f, scalar &ePore, scalar &EigenError) const;

        //- Calculate Eigen values and vectors
        void calculateEigens(
            vector &sigma_prin,
            tensor &ev,
            const symmTensor sigma,
            scalar &EigenError) const;

        //- Disallow default bitwise copy construct
        linearElasticMohrCoulombPlasticDilationCutoff(const linearElasticMohrCoulombPlasticDilationCutoff &);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticMohrCoulombPlasticDilationCutoff &);

    public:
        //- Runtime type information
        TypeName("linearElasticMohrCoulombPlasticDilationCutoff");

        // Static data members

        // Constructors

        //- Construct from dictionary
        linearElasticMohrCoulombPlasticDilationCutoff(
            const word &name,
            const fvMesh &mesh,
            const dictionary &dict,
            const nonLinearGeometry::nonLinearType &nonLinGeom);

        // Destructor

        virtual ~linearElasticMohrCoulombPlasticDilationCutoff();

        // Member Functions

        //- Return density
        //(virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar &mu() const;

        //- Provide access to lambda for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar &lambda() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField &sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField &sigma);

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
