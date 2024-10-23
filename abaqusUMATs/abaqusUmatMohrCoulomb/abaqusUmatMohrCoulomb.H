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
    abaqusUmatMohrCoulomb

Description
    Wrapper class for abaqusUmatMohrCoulomb.f fortran sub-routine from Abaqus.

SourceFiles
    abaqusUmatMohrCoulomb.C

Author
    Denis Maier, BAW.

    Based on the work by Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef abaqusUmatMohrCoulomb_H
#define abaqusUmatMohrCoulomb_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                         Class abaqusUmatMohrCoulomb Declaration
\*---------------------------------------------------------------------------*/

    class abaqusUmatMohrCoulomb
        : public mechanicalLaw
    {
        // Private data

        //- Density --> now in materialLaw
        //const dimensionedScalar rho_;

        //- List of material properties
        double PROPS[5];

        dimensionedScalar E_;
        dimensionedScalar nu_;
        dimensionedScalar lambda_;
        dimensionedScalar mu_;
        dimensionedScalar K_;
        dimensionedScalar varPhi_;
        dimensionedScalar c_;
        dimensionedScalar varPsi_;
        //- List of internal state scalar variable fields
        // PtrList<volScalarField> stateVariables_;

        //- "Implicit stiffness" used by the segregated solid models
        //  This only affects the convergence assuming convergence is achieved
        //  For linear elastic solids the ideal value is "2*mu + lambda"
        const dimensionedScalar impK_;

        volScalarField DepsilonDevEq_;

        volScalarField activeYield_;

        scalar resSig_;

        // Private Member Functions

        //- Disallow default bitwise copy construct
        abaqusUmatMohrCoulomb(const abaqusUmatMohrCoulomb &);

        //- Disallow default bitwise assignment
        void operator=(const abaqusUmatMohrCoulomb &);

    public:
        //- Runtime type information
        TypeName("abaqusUmatMohrCoulomb");

        // Constructors

        //- Construct from dictionary
        abaqusUmatMohrCoulomb(
            const word &name,
            const fvMesh &mesh,
            const dictionary &dict,
            const nonLinearGeometry::nonLinearType &nonLinGeom);

        // Destructor

        virtual ~abaqusUmatMohrCoulomb();

        // Member Functions

        //- Return density --> now in materialLaw
        //virtual tmp<volScalarField> rho() const;

        //- Return density scale --> now in materialLaw
        //virtual const dimensionedScalar &rhoScalar() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar &mu() const;

        //- Lambda
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
