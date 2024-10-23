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
    MFront

Description
    Wrapper class for MFront. UNFINISHED!

SourceFiles
    MFront.C

Author
    Denis Maier, BAW. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MFront_H
#define MFront_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                         Class MFront Declaration
\*---------------------------------------------------------------------------*/

    class MFront
        : public mechanicalLaw
    {
        // Private data
        const word MFrontLibPath;
        const word behaviourName;
        mgis::behaviour::Hypothesis hypothesis_;
        mgis::behaviour::Behaviour behaviour_;
        //int const equivalent_plastic_strain_offset_;

        //- List of material properties
        dictionary properties_;

        dimensionedScalar rho_;
        //- Total strain
        volSymmTensorField epsilon_;

        mgis::behaviour::MaterialDataManager MFrontBehaviourData;
        // Private Member Functions

        //- Disallow default bitwise copy construct
        MFront(const MFront &);

        //- Disallow default bitwise assignment
        void operator=(const MFront &);

        const word *toString(mgis::behaviour::Behaviour::Kinematic);
        const word *toString(mgis::behaviour::Behaviour::Symmetry);
        const word *bTypeToString(int);
        const word *varTypeToString(int);

    public:
        //- Runtime type information
        TypeName("MFront");

        // Constructors

        //- Construct from dictionary
        MFront(
            const word &name,
            const fvMesh &mesh,
            const dictionary &dict,
            const nonLinearGeometry::nonLinearType &nonLinGeom);

        // Destructor

        virtual ~MFront();

        // Member Functions

        void updateParams();

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return density scale
        virtual const dimensionedScalar &rhoScalar() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar &mu() const;

        //- Bulk modulus
        virtual const dimensionedScalar &K() const;

        //- Young's modulus
        virtual const dimensionedScalar &E() const;

        //- Poisson's ratio
        virtual const dimensionedScalar &nu() const;

        //- Lambda
        virtual const dimensionedScalar &lambda() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField &sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField &sigma);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
