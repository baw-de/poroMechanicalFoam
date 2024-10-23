/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fixedPoroPotentialFvPatchScalarField

Description
    Foam::fixedPoroPotentialFvPatchScalarField

SourceFiles
    fixedPoroPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPoroPotentialFvPatchScalarField_H
#define fixedPoroPotentialFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "PatchFunction1.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class totalPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedPoroPotentialFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Total head
        scalarField h0_;

        //- head time series
        autoPtr<PatchFunction1<scalar>> headSeries_;

        //- is field pressure head or pressure
        Switch isHead_;

        //- where to look for fields when hydro-mechanically coupled
        word HMCoupled;

public:

    //- Runtime type information
    TypeName("fixedPoroPotential");


    // Constructors

        //- Construct from patch and internal field
        fixedPoroPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPoroPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedPoroPotentialFvPatchScalarField
        //  onto a new patch
        fixedPoroPotentialFvPatchScalarField
        (
            const fixedPoroPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedPoroPotentialFvPatchScalarField
        (
            const fixedPoroPotentialFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedPoroPotentialFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedPoroPotentialFvPatchScalarField
        (
            const fixedPoroPotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedPoroPotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the total head
            const scalarField& h0() const
            {
                return h0_;
            }

            //- Return reference to the total head to allow adjustment
            scalarField& h0()
            {
                return h0_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
