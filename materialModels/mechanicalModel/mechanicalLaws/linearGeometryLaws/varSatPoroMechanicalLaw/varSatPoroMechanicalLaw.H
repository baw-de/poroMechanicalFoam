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
    varSatPoroMechanicalLaw

Description
    Generalised poro-coupled mechanical law, where the effective stress is
    defined by any other run-time selectable mechanical law, to which a pore
    stress contribution is added.

    An optional Biot's coefficient can be given to scale the pressure term.

    The effect of partial saturation is added through effective stress models,
    scaling the pressure contribution to effective stress by some function of
    the degree of saturation.

SourceFiles
    varSatPoroMechanicalLaw.C

Author
    Denis Maier, BAW. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroMechanicalLaw_H
#define varSatPoroMechanicalLaw_H

#include "mechanicalLaw.H"
#include "scalarDiskReader.H"
#include "effectiveStressModel.H"
#include "poroMechanicalLaw2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class varSatPoroMechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class varSatPoroMechanicalLaw
:
    public poroMechanicalLaw2
{
    // Private Data

        //- scalar reader function to read saturation from file
        mutable autoPtr<scalarDiskReader> saturationReaderPtr_;

        //- Run-time selectable mechanical law to define the effective stress
        autoPtr<effectiveStressModel> effectiveStressModelPtr_;

        //- If pressure and saturation were read from disk, we need to read porosity
        // from dict.
        mutable autoPtr<const volScalarField> nSubMesh_;

        //- If pressure and saturation were read from disk, we need to read water density
        // from dict.
        const dimensionedScalar rho_water;

    // Private Member Functions

        scalarDiskReader &saturationReader() const
        {
            return saturationReaderPtr_.ref();
        }

        //- check if sigmaEff has been initalized already
        bool checkSigmaEffReady(const volSymmTensorField& sigma, const volScalarField& p, const volScalarField& chi);

        //- check if sigmaEfff has been initalized already
        bool checkSigmaEffReady(const surfaceSymmTensorField& sigma, const surfaceScalarField& p, const surfaceScalarField& chi);

        //- return reference to the effectiveStressModel / effective stress definition
        effectiveStressModel &effectiveStressModelRef();

        tmp<volScalarField> n() const
        {
            if(!nSubMesh_.valid())
            {
                makeN();
            }
            return nSubMesh_();
        }

        void makeN() const;

        //- Disallow default bitwise copy construct
        varSatPoroMechanicalLaw(const varSatPoroMechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const varSatPoroMechanicalLaw&);

public:

    //- Runtime type information
    TypeName("varSatPoroMechanicalLaw");

    // Static data members


    // Constructors

        //- Construct from dictionary
        varSatPoroMechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~varSatPoroMechanicalLaw();


    // Member Functions

        //- Return updated total density
        virtual tmp<volScalarField> rho() const;

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
