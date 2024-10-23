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

#include "zDependendDisplacementOrTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zDependendDisplacementOrTractionFvPatchVectorField::
zDependendDisplacementOrTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    constantDisplacement_(p.size(), vector::zero),
    constantTraction_(p.size(), vector::zero),
    constantPressure_(p.size(), 0.0),
    z0_(0.0),
    splitDir_(0),
    displacementSeries_(),
    tractionSeries_(),
    pressureSeries_(),
    zSeries_()
{}


zDependendDisplacementOrTractionFvPatchVectorField::
zDependendDisplacementOrTractionFvPatchVectorField
(
    const zDependendDisplacementOrTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    constantDisplacement_(ptf.constantDisplacement_, mapper),
    constantTraction_(ptf.constantTraction_, mapper),
    constantPressure_(ptf.constantPressure_, mapper),
    z0_(ptf.z0_),
    splitDir_(ptf.splitDir_),
    displacementSeries_(ptf.displacementSeries_),
    tractionSeries_(ptf.tractionSeries_),
    pressureSeries_(ptf.pressureSeries_),
    zSeries_(ptf.zSeries_)
{}


zDependendDisplacementOrTractionFvPatchVectorField::
zDependendDisplacementOrTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    constantDisplacement_(p.size(), vector::zero),
    constantTraction_(p.size(), vector::zero),
    constantPressure_(p.size(), 0.0),
    z0_(0.0),
    splitDir_(0),
    displacementSeries_(),
    tractionSeries_(),
    pressureSeries_(),
    zSeries_()
{

    splitDir_ = Switch(dict.lookup("splitDirection"));

    // Check if displacement is time-varying
    if (dict.found("displacementSeries") && dict.found("constantDisplacement"))
    {
        FatalErrorIn
        (
            "zDependendDisplacementOrTractionFvPatchVectorField::"
            "zDependendDisplacementOrTractionFvPatchVectorField"
        )   << "constantDisplacement or displacementSeries can be specified, "
            << "not both!" << abort(FatalError);
    }
    else if (dict.found("displacementSeries"))
    {
        Info<< type() << ": " << patch().name()
            << " displacement is time-varying" << endl;
        displacementSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));
    }
    else
    {
        constantDisplacement_ = vectorField("constantDisplacement", dict, p.size());
    }

    // Check if traction is time-varying
    if (dict.found("tractionSeries") && dict.found("constantTraction"))
    {
        FatalErrorIn
        (
            "zDependendDisplacementOrTractionFvPatchVectorField::"
            "zDependendDisplacementOrTractionFvPatchVectorField"
        )   << "constantTraction or tractionSeries can be specified, "
            << "not both!" << abort(FatalError);
    }
    else if (dict.found("tractionSeries"))
    {
        Info<< type() << ": " << patch().name()
            << " traction is time-varying" << endl;
            tractionSeries_ =
                interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else
    {
        constantTraction_ = vectorField("constantTraction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries") && dict.found("constantPressure"))
    {
        FatalErrorIn
        (
            "zDependendDisplacementOrTractionFvPatchVectorField::"
            "zDependendDisplacementOrTractionFvPatchVectorField"
        )   << "constantPressure or pressureSeries can be specified, "
            << "not both!" << abort(FatalError);
    }
    else if (dict.found("pressureSeries"))
    {
        Info<< type() << ": " << patch().name()
            << " pressure is time-varying" << endl;
    }
    else
    {
        constantPressure_ = scalarField("constantPressure", dict, p.size());
    }


            // Check if z0_ is time-varying
            if (dict.found("zSeries") && dict.found("z0"))
            {
                FatalErrorIn
                (
                    "zDependendDisplacementOrTractionFvPatchVectorField::"
                    "zDependendDisplacementOrTractionFvPatchVectorField"
                )   << "z0 or zSeries can be specified, "
                    << "not both!" << abort(FatalError);
            }
            else if (dict.found("zSeries"))
            {
                Info<< type() << ": " << patch().name()
                    << " z0 is time-varying" << endl;
                    zSeries_ =
                        interpolationTable<scalar>(dict.subDict("zSeries"));
            }
            else
            {
                z0_ = readScalar(dict.lookup("z0"));
            }




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
        FatalErrorIn
        (
            "zDependendDisplacementOrTractionFvPatchVectorField::"
            "zDependendDisplacementOrTractionFvPatchVectorField"
        )   << "value or refValue entry must be specified for patch "
            << patch().name() << abort(FatalError);
    }

    if (zSeries_.size())
    {
      z0_ = zSeries_(db().time().timeOutputValue());
    }

    valueFraction() = sqr(patch().nf());
    refGrad() = vectorField(patch().size(), vector::zero);
    refValue() = constantDisplacement_;

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        const Field<vector> normalValue =
            transform(valueFraction(), refValue());

        const Field<vector> gradValue =
            patchInternalField() + refGrad()/patch().deltaCoeffs();

        const Field<vector> transformGradValue =
            transform(I - valueFraction(), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }

}


zDependendDisplacementOrTractionFvPatchVectorField::
zDependendDisplacementOrTractionFvPatchVectorField
(
    const zDependendDisplacementOrTractionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    constantDisplacement_(ptf.constantDisplacement_),
    constantTraction_(ptf.constantTraction_),
    constantPressure_(ptf.constantPressure_),
    z0_(ptf.z0_),
    splitDir_(ptf.splitDir_),
    displacementSeries_(ptf.displacementSeries_),
    tractionSeries_(ptf.tractionSeries_),
    pressureSeries_(ptf.pressureSeries_),
    zSeries_(ptf.zSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void zDependendDisplacementOrTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);

    constantDisplacement_.autoMap(m);
    constantTraction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void zDependendDisplacementOrTractionFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const zDependendDisplacementOrTractionFvPatchVectorField& dmptf =
        refCast<const zDependendDisplacementOrTractionFvPatchVectorField>(ptf);

    constantDisplacement_.rmap(dmptf.constantDisplacement_, addr);
    constantTraction_.rmap(dmptf.constantTraction_, addr);
}


void zDependendDisplacementOrTractionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Note: valueFraction is set in the constructor so no need to set it again
  if (zSeries_.size())
  {
    z0_ = zSeries_(db().time().timeOutputValue());
  }

  const fvPatchField<scalar>& z_ =
      patch().lookupPatchField<volScalarField, scalar>("z");

      const vectorField nf = patch().nf();

  forAll(valueFraction(),faceI)
  {
    if (z_[faceI] < z0_)
    {
      valueFraction()[faceI] = sqr(nf[faceI]) * splitDir_;
    }
    else
    {
      valueFraction()[faceI] = sqr(nf[faceI]) * (1 - splitDir_);
    }
  }



    // Get current displacement field to be specified
    vectorField disp = constantDisplacement_;

    if (displacementSeries_.size())
    {
        disp = displacementSeries_(db().time().timeOutputValue());
    }

    // Lookup the solidModel object
    const solidModel& solMod =
        lookupSolidModel(patch().boundaryMesh().mesh());

    if (solMod.incremental())
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    // Set displacement
    refValue() = disp;

    // Unit normals
    const vectorField Sf = patch().Sf();

    // Get current traction to be specified (defaults to zero)
    vectorField traction = constantTraction_ + Sf * constantPressure_;

if (pressureSeries_.size() && tractionSeries_.size())
{
  traction = tractionSeries_(db().time().timeOutputValue())
              - Sf * pressureSeries_(db().time().timeOutputValue());
}
else if (pressureSeries_.size() || tractionSeries_.size())
{
    if (tractionSeries_.size())
    {
        traction = tractionSeries_(db().time().timeOutputValue()) - Sf * constantPressure_ ;
    }

    if (pressureSeries_.size())
    {
        traction = constantTraction_ - Sf * pressureSeries_(db().time().timeOutputValue());
    }
}
    // Iteratively set gradient to enforce the specified traction
    refGrad() =
        solMod.tractionBoundarySnGrad
        (
            traction,
            scalarField(patch().size(), 0.0),
            patch()
        );

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void zDependendDisplacementOrTractionFvPatchVectorField::write(Ostream& os) const
{
    if (displacementSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        displacementSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        constantDisplacement_.writeEntry("constantDisplacement", os);
    }

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        constantTraction_.writeEntry("constantTraction", os);
    }


    if (pressureSeries_.size())
    {
        os.writeKeyword("presureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        constantPressure_.writeEntry("constantPressure", os);
    }

    if (zSeries_.size())
    {
        os.writeKeyword("zSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        zSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        os << tab << "z0" << tab << z0_ << ";" << endl;
    }

    os << tab << "splitDirection" << tab << splitDir_ << ";" << endl;

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    zDependendDisplacementOrTractionFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
