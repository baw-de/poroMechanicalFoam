/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    poroMechanicalLaw2

Description
    poroMechanicalLaw from solids4Foam but
    - added the ability to read pressure from disk
    - works more conveniently with poroMechanicalFoam

SourceFiles
    poroMechanicalLaw2.C

Authors
    Heavily leaning on code from Philip Cardiff, UCD. All rights reserved.
    Denis Maier, BAW.

\*---------------------------------------------------------------------------*/

#ifndef poroMechanicalLaw_H
#define poroMechanicalLaw_H

#include "mechanicalLaw.H"
#include "scalarDiskReader.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poroMechanicalLaw2 Declaration
\*---------------------------------------------------------------------------*/

class poroMechanicalLaw2
:
    public mechanicalLaw
{
    // Private Data
        protected:

        //- Run-time selectable mechanical law to define the effective stress
        autoPtr<mechanicalLaw> effectiveStressMechLawPtr_;

        //- effective stress
        autoPtr<volSymmTensorField> sigmaEff_;

        //- effective stress on surfaceMesh
        autoPtr<surfaceSymmTensorField> sigmaEfff_;

        const Switch buoyancySwitch_;

        //- Name of pressure field (defaults to "p")
        const word pName_;

        const Switch readFromDisk_;

        autoPtr<scalarDiskReader> pReaderPtr_;

        //- Biot's coefficient
        //  Optional coefficient to scalar the pressure contribution to the
        //  total stress. Defaults to 1.0.
        const dimensionedScalar b_;

    // Private Member Functions

        //- check if sigmaEff has been initalized already
        bool checkSigmaEffReady(const volSymmTensorField& sigma, const volScalarField& p);

        //- check if sigmaEfff has been initalized already
        bool checkSigmaEffReady(const surfaceSymmTensorField& sigma, const surfaceScalarField& p);

        //- Make the initial pore-pressure surface field
        void makeP0f() const;

        //- Return a reference to the initial pore-pressure surface field
        const surfaceScalarField& p0f() const;

        //- Disallow default bitwise copy construct
        poroMechanicalLaw2(const poroMechanicalLaw2&);

        //- Disallow default bitwise assignment
        void operator=(const poroMechanicalLaw2&);

protected:
        //- Lookup and return a reference to the pressure field
        tmp<volScalarField> lookupFluidField(const word fieldName) const;

        scalarDiskReader& pReader()
        {
            return pReaderPtr_.ref();
        }

public:

    //- Runtime type information
    TypeName("poroMechanicalLaw2");

    // Static data members


    // Constructors

        //- Construct from dictionary
        poroMechanicalLaw2
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~poroMechanicalLaw2();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        virtual tmp<volScalarField> biotCoeff() const
        {
            tmp<volScalarField> tb(
                new volScalarField(
                  IOobject(
                  "biotCoeff",
                  mesh().time().timeName(),
                  mesh(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
                  ),
                  mesh(),
                  b_
                )
            );
            return tb;
        }

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
