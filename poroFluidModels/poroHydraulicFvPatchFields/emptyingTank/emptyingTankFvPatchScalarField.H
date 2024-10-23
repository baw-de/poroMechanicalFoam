/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::emptyingTankFvPatchScalarField

Description
    Explicit boundary condition for tank of fluid infiltrating in porous Media and appropriate head reduction.

SourceFiles
    emptyingTankFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef emptyingTankFvPatchScalarField_H
#define emptyingTankFvPatchScalarField_H

#ifdef FOAMEXTEND
#include "foamTime.H"
#endif
#include "mixedFvPatchFields.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
           Class emptyingTankFvPatchField Declaration
\*---------------------------------------------------------------------------*/

    class emptyingTankFvPatchScalarField
        : public mixedFvPatchScalarField
    {
        // Private data

        //- Name of the z field
        word pName_;

        //- Name of the flux transporting the field
        word fluxname_;

        //- Total temperature
        scalar h0_;
        //- Total temperature
        scalar dh_;
        //- Pressure time series
        dimensionedScalar time_;
        //- Pressure time series
        interpolationTable<scalar> crossSectionSeries_;
        //- where to look for fields when hydro-mechanically coupled
        word HMCoupled;

    public:
        //- Runtime type information
        TypeName("emptyingTank");

        // Constructors

        //- Construct from patch and internal field
        emptyingTankFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        emptyingTankFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  emptyingTankFvPatchScalarField onto a new patch
        emptyingTankFvPatchScalarField(
            const emptyingTankFvPatchScalarField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        emptyingTankFvPatchScalarField(
            const emptyingTankFvPatchScalarField &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new emptyingTankFvPatchScalarField(*this));
        }

        //- Construct as copy setting internal field reference
        emptyingTankFvPatchScalarField(
            const emptyingTankFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new emptyingTankFvPatchScalarField(*this, iF));
        }

        // Member functions

        // Access

        //- Return the total pressure
        const scalar &h0() const
        {
            return h0_;
        }

        //- Return reference to the total pressure to allow adjustment
        scalar &h0()
        {
            return h0_;
        }

        // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(
            const fvPatchFieldMapper &);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap(
            const fvPatchScalarField &,
            const labelList &);

        // Evaluation functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream &) const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
