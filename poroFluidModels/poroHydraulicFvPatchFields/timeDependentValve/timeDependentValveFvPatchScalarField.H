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
    Foam::timeDependentValveFvPatchScalarField

Description
    Foam::timeDependentValveFvPatchScalarField

SourceFiles
    timeDependentValveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeDependentValveFvPatchScalarField_H
#define timeDependentValveFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "interpolationTable.H"
#include "patchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
           Class timeDependentOutletWallFvPatchField Declaration
\*---------------------------------------------------------------------------*/

    class timeDependentValveFvPatchScalarField
        : public mixedFvPatchScalarField
    {
        // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Heat capacity ratio
        scalarField outletValue_;

        //- outlet valve time series
        autoPtr<patchFunction1> valveStateSeries;

    public:
        //- Runtime type information
        TypeName("timeDependentOutletWall");

        // Constructors

        //- Construct from patch and internal field
        timeDependentValveFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        timeDependentValveFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  timeDependentValveFvPatchScalarField onto a new patch
        timeDependentValveFvPatchScalarField(
            const timeDependentValveFvPatchScalarField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        timeDependentValveFvPatchScalarField(
            const timeDependentValveFvPatchScalarField &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new timeDependentValveFvPatchScalarField(*this));
        }

        //- Construct as copy setting internal field reference
        timeDependentValveFvPatchScalarField(
            const timeDependentValveFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new timeDependentValveFvPatchScalarField(*this, iF));
        }

        // Member functions

        // Access

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
