/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::fixedPoroFluxFvPatchScalarField

Description
    Foam::fixedPoroFluxFvPatchScalarField

SourceFiles
    fixedPoroFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPoroFluxFvPatchScalarFields_H
#define fixedPoroFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fixedPoroFluxFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedPoroFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data
    scalarField flux_;

    //- flux time series
    autoPtr<PatchFunction1<scalar>> fluxSeries_;

    //- is the field pressure head or pressure?
    Switch isHead_;

    //- where to look for fields when hydro-mechanically coupled
    word HMCoupled;

public:

    //- Runtime type information
    TypeName("fixedPoroFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedPoroFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPoroFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedPoroFluxFvPatchScalarField onto
        //  a new patch
        fixedPoroFluxFvPatchScalarField
        (
            const fixedPoroFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedPoroFluxFvPatchScalarField
        (
            const fixedPoroFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedPoroFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedPoroFluxFvPatchScalarField
        (
            const fixedPoroFluxFvPatchScalarField&,
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
                new fixedPoroFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
