/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "fixedWallNoTension.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    fixedWallNoTension::
        fixedWallNoTension(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF)
        : solidDirectionMixedFvPatchVectorField(p, iF)
    {
    }

    fixedWallNoTension::
        fixedWallNoTension(
            const fixedWallNoTension &ptf,
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const fvPatchFieldMapper &mapper)
        : solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper)
    {
    }

    fixedWallNoTension::
        fixedWallNoTension(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const dictionary &dict)
        : solidDirectionMixedFvPatchVectorField(p, iF)
    {

        // Lookup and set refValue, refGrad, valueFraction and value
        if (dict.found("refValue"))
        {
            refValue() = vectorField("refValue", dict, p.size());
        }
        else if (dict.found("value"))
        {
            refValue() = vectorField("value", dict, p.size());
        }
        else
        {
            FatalErrorIn(
                "fixedWallNoTension::"
                "fixedWallNoTension")
                << "value or refValue entry must be specified for patch "
                << patch().name() << abort(FatalError);
        }

        valueFraction() = sqr(patch().nf());
        refGrad() = vectorField(patch().size(), vector::zero);
        refValue() = vectorField(patch().size(), vector::zero);

        if (dict.found("value"))
        {
            Field<vector>::operator=(vectorField("value", dict, p.size()));
        }
        else
        {
            const Field<vector> normalValue =
                transform(valueFraction(), refValue());

            const Field<vector> gradValue =
                patchInternalField() + refGrad() / patch().deltaCoeffs();

            const Field<vector> transformGradValue =
                transform(I - valueFraction(), gradValue);

            Field<vector>::operator=(normalValue + transformGradValue);
        }
    }

    fixedWallNoTension::
        fixedWallNoTension(
            const fixedWallNoTension &ptf,
            const DimensionedField<vector, volMesh> &iF)
        : solidDirectionMixedFvPatchVectorField(ptf, iF)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    // Map from self
    void fixedWallNoTension::autoMap(
        const fvPatchFieldMapper &m)
    {
        solidDirectionMixedFvPatchVectorField::autoMap(m);
    }

    // Reverse-map the given fvPatchField onto this fvPatchField
    void fixedWallNoTension::rmap(
        const fvPatchField<vector> &ptf,
        const labelList &addr)
    {
        solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);
        const fixedWallNoTension &dmptf =
            refCast<const fixedWallNoTension>(ptf);
    }

    void fixedWallNoTension::updateCoeffs()
    {
        if (this->updated())
        {
            return;
        }

        // Note: valueFraction is set in the constructor so no need to set it again

        const fvPatchField<vector> &D_ =
            patch().lookupPatchField<volVectorField, vector>("D");

        const vectorField nf = patch().nf();

        scalarField a = (D_ / max(SMALL, mag(D_))) & nf;

        Info << "min a is: " << min(a) << endl;
        forAll(valueFraction(), faceI)
        {
            if (a[faceI] > 0)
            {
                valueFraction()[faceI] = sqr(nf[faceI]);
            }
            else
            {
                valueFraction()[faceI] = sqr(nf[faceI]) * 0;
            }
        }

        // Lookup the solidModel object
        const solidModel &solMod =
            lookupSolidModel(patch().boundaryMesh().mesh());

        // Get current traction to be specified (defaults to zero)
        vectorField traction(patch().size(), vector::zero);

        // Iteratively set gradient to enforce the specified traction
        refGrad() =
            solMod.tractionBoundarySnGrad(
                traction,
                scalarField(patch().size(), 0.0),
                patch());

        solidDirectionMixedFvPatchVectorField::updateCoeffs();
    }

    // Write
    void fixedWallNoTension::write(Ostream &os) const
    {
        solidDirectionMixedFvPatchVectorField::write(os);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField(
        fvPatchVectorField,
        fixedWallNoTension);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
