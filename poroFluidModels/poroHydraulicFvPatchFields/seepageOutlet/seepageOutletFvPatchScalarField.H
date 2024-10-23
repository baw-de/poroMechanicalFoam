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
    Foam::seepageOutletFvPatchScalarField

Description
    Foam::seepageOutletFvPatchScalarField

SourceFiles
    seepageOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef seepageOutletFvPatchScalarField_H
#define seepageOutletFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "PatchFunction1.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
           Class seepageOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

    class seepageOutletFvPatchScalarField
        : public mixedFvPatchScalarField
    {
        // Private data

        //- Name of the z field
        word zName_;

        //- Heat capacity ratio
        scalarField outletPressure_;

        //- Total temperature
        scalarField h0_;

        //- Pressure time series
        autoPtr<PatchFunction1<scalar>> hSeries_;

        // Use the Depuit approximation?
        Switch depuitApprox_;

        //- is field pressure head or pressure
        Switch isHead_;

        //- where to look for fields when hydro-mechanically coupled
        word HMCoupled;

        //- update this bc explicitly after each timestep?
        Switch explicit_;

        //- time the boundary was updated the last time
        scalar updatedTime_;

    public:
        //- Runtime type information
        TypeName("seepageOutlet");

        // Constructors

        //- Construct from patch and internal field
        seepageOutletFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        seepageOutletFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  seepageOutletFvPatchScalarField onto a new patch
        seepageOutletFvPatchScalarField(
            const seepageOutletFvPatchScalarField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        seepageOutletFvPatchScalarField(
            const seepageOutletFvPatchScalarField &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new seepageOutletFvPatchScalarField(*this));
        }

        //- Construct as copy setting internal field reference
        seepageOutletFvPatchScalarField(
            const seepageOutletFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new seepageOutletFvPatchScalarField(*this, iF));
        }

        // Member functions

        // Access

        //- Return the total pressure
        const scalarField &h0() const
        {
            return h0_;
        }

        //- Return reference to the total pressure to allow adjustment
        scalarField &h0()
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
