/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    ohdeElastic

Description
    Hookean linear elastic mechanical law.

SourceFiles
    ohdeElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef ohdeElastic_H
#define ohdeElastic_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ohdeElastic Declaration
\*---------------------------------------------------------------------------*/

class ohdeElastic
:
    public mechanicalLaw
{
    // Private data

        //- Shear modulus (second Lame parameter)
        dimensionedScalar mu0_;

        //- Bulk modulus
        dimensionedScalar K0_;
        
        dimensionedScalar n_;
        
        dimensionedScalar p0_;

        //- Young's modulus
        volScalarField K_;

        //- Poisson's ratio
        volScalarField mu_;
        
        //- Young's modulus
        volScalarField E_;

        //- Poisson's ratio
        volScalarField nu_;

        //- First Lame parameter
        volScalarField lambda_;


    // Private Member Functions

        void updateStiffness(const volSymmTensorField& sigma);
        
        //- Disallow default bitwise copy construct
        ohdeElastic(const ohdeElastic&);

        //- Disallow default bitwise assignment
        void operator=(const ohdeElastic&);

public:

    //- Runtime type information
    TypeName("ohdeElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        ohdeElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~ohdeElastic()
        {}


    // Member Functions

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress given epsilon
        virtual void correct
        (
            volSymmTensorField& sigma, const volSymmTensorField& epsilon
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
