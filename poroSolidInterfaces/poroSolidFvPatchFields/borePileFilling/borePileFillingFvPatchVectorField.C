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

#include "borePileFillingFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

borePileFillingFvPatchVectorField::
borePileFillingFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    constantDisplacement_(p.size(),vector::zero),
    z0_(0.0),
    h_c_(0.0),
    gamma_c_(25000),
    splitDir_(0),
    hcSeries_(),
    zSeries_()
{}


borePileFillingFvPatchVectorField::
borePileFillingFvPatchVectorField
(
    const borePileFillingFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    constantDisplacement_(ptf.constantDisplacement_),
    z0_(ptf.z0_),
    h_c_(ptf.h_c_),
    gamma_c_(ptf.gamma_c_),
    splitDir_(ptf.splitDir_),
    hcSeries_(ptf.hcSeries_),
    zSeries_(ptf.zSeries_)
{}


borePileFillingFvPatchVectorField::
borePileFillingFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    constantDisplacement_(p.size(),vector::zero),
    z0_(0.0),
    h_c_(0.0),
    gamma_c_(25000),
    splitDir_(0),
    hcSeries_(),
    zSeries_()
{

    splitDir_ = Switch(dict.lookup("splitDirection"));


            // Check if z0_ is time-varying
            if (dict.found("zSeries") && dict.found("z0"))
            {
                FatalErrorIn
                (
                    "borePileFillingFvPatchVectorField::"
                    "borePileFillingFvPatchVectorField"
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

            // Check if h_c_ is time-varying
            if (dict.found("hcSeries") && dict.found("hc"))
            {
                FatalErrorIn
                (
                    "borePileFillingFvPatchVectorField::"
                    "borePileFillingFvPatchVectorField"
                )   << "hc or hcSeries can be specified, "
                    << "not both!" << abort(FatalError);
            }
            else if (dict.found("hcSeries"))
            {
                Info<< type() << ": " << patch().name()
                    << " hc is time-varying" << endl;
                    hcSeries_ =
                        interpolationTable<scalar>(dict.subDict("hcSeries"));
            }
            else
            {
                h_c_ = readScalar(dict.lookup("hc"));
            }

            if (dict.found("gamma_c"))
            {
                gamma_c_ = readScalar(dict.lookup("gamma_c"));
            }

            Info << "Specific Weight of concrete is: "  << gamma_c_ << endl;

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
            "borePileFillingFvPatchVectorField::"
            "borePileFillingFvPatchVectorField"
        )   << "value or refValue entry must be specified for patch "
            << patch().name() << abort(FatalError);
    }

    if (zSeries_.size())
    {
      z0_ = zSeries_(db().time().timeOutputValue());
    }
    if (hcSeries_.size())
    {
      h_c_ = hcSeries_(db().time().timeOutputValue());
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


borePileFillingFvPatchVectorField::
borePileFillingFvPatchVectorField
(
    const borePileFillingFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    constantDisplacement_(ptf.constantDisplacement_),
    z0_(ptf.z0_),
    h_c_(ptf.h_c_),
    gamma_c_(ptf.gamma_c_),
    splitDir_(ptf.splitDir_),
    hcSeries_(ptf.hcSeries_),
    zSeries_(ptf.zSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void borePileFillingFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);

    constantDisplacement_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void borePileFillingFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const borePileFillingFvPatchVectorField& dmptf =
        refCast<const borePileFillingFvPatchVectorField>(ptf);

    constantDisplacement_.rmap(dmptf.constantDisplacement_, addr);
}


void borePileFillingFvPatchVectorField::updateCoeffs()
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

  if (hcSeries_.size())
  {
    h_c_ = hcSeries_(db().time().timeOutputValue());
  }

  const fvPatchField<scalar>& z_ =
      patch().lookupPatchField<volScalarField, scalar>("z");

      // Unit normals
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
    refValue()  = constantDisplacement_;

    // Lookup the solidModel object
    const solidModel& solMod =
        lookupSolidModel(patch().boundaryMesh().mesh());

    const fvPatchField<scalar>& p_ =
        patch().lookupPatchField<volScalarField, scalar>("pHead");

    const uniformDimensionedVectorField& gamma_w_ =
            db().lookupObject<uniformDimensionedVectorField>("gamma_w");

    // Get current traction to be specified (defaults to zero)
    scalarField effPressure = (h_c_-z_)*gamma_c_ - p_ * mag(gamma_w_.value());

    Info << "max traction on boundary is: " << max(mag(effPressure)) << endl
         << "with effective stress: " << max(mag(nf * ((h_c_-z_)*gamma_c_))) << endl
         << "and pore pressure: " << max(mag(nf * (p_ * mag(gamma_w_.value())))) << endl;
    // Iteratively set gradient to enforce the specified traction
    refGrad() =
        solMod.tractionBoundarySnGrad
        (
            vectorField(patch().size(),vector::zero),
            effPressure,
            patch()
        );

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void borePileFillingFvPatchVectorField::write(Ostream& os) const
{

        constantDisplacement_.writeEntry("constantDisplacement", os);

    if (hcSeries_.size())
    {
        os.writeKeyword("hcSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        hcSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
      os << tab << "hc" << tab << h_c_ << ";" << endl;
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
    borePileFillingFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
