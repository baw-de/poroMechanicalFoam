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

\*---------------------------------------------------------------------------*/

#include "timeDependentValveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeDependentValveFvPatchScalarField::
    timeDependentValveFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(p, iF),
      phiName_("phi"),
      outletValue_(0.0),
      valveStateSeries()
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}

Foam::timeDependentValveFvPatchScalarField::
    timeDependentValveFvPatchScalarField(
        const timeDependentValveFvPatchScalarField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : mixedFvPatchScalarField(ptf, p, iF, mapper),
      phiName_(ptf.phiName_),
      outletValue_(ptf.outletValue_),
      valveStateSeries(ptf.valveStateSeries.clone(this->patch().patch()))
{
}

Foam::timeDependentValveFvPatchScalarField::
    timeDependentValveFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
    : mixedFvPatchScalarField(p, iF),
      phiName_(dict.lookupOrDefault<word>("flux", "phi")),
      outletValue_("outletValue", dict, p.size()),
      valveStateSeries(PatchFunction1<scalar>::New(p.patch(), "valveStateSeries", dict))
{
    this->refValue() = pTraits<scalar>::zero;
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=(
            scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchField<scalar>::operator=(outletValue_);
    }

    this->refGrad() = scalarField(this->size(), 0.0);
    this->valueFraction() = 0.0;
}

Foam::timeDependentValveFvPatchScalarField::
    timeDependentValveFvPatchScalarField(
        const timeDependentValveFvPatchScalarField &tppsf)
    : mixedFvPatchScalarField(tppsf),
      phiName_(tppsf.phiName_),
      outletValue_(tppsf.outletValue_),
      valveStateSeries(tppsf.valveStateSeries.clone(this->patch().patch()))
{
}

Foam::timeDependentValveFvPatchScalarField::
    timeDependentValveFvPatchScalarField(
        const timeDependentValveFvPatchScalarField &tppsf,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(tppsf, iF),
      phiName_(tppsf.phiName_),
      outletValue_(tppsf.outletValue_),
      valveStateSeries(tppsf.valveStateSeries.clone(this->patch().patch()))
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeDependentValveFvPatchScalarField::autoMap(
    const fvPatchFieldMapper &m)
{
    mixedFvPatchScalarField::autoMap(m);
}

void Foam::timeDependentValveFvPatchScalarField::rmap(
    const fvPatchScalarField &ptf,
    const labelList &addr)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const timeDependentValveFvPatchScalarField &tiptf =
        refCast<const timeDependentValveFvPatchScalarField>(ptf);
}

void Foam::timeDependentValveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar valveState = valveStateSeries(this->db().time().timeOutputValue());

    const poroHydraulicModel *poroHydraulic_ = &this->db().time().lookupObject<poroHydraulicModel>("poroHydraulicModel");
    scalarField n_ = patch().nf() & vector(poroHydraulic_->gamma().value()).normalise();

    scalarField tmpp(outletValue_);
    scalarField tmpValFrac(this->size(), 0.0);

    scalarField pInternal = patchInternalField();

    if (valveState < 1.0)
    {
        forAll(tmpp, iface)
        {
            tmpp[iface] = pInternal[iface];
            tmpValFrac[iface] = 0.0;
        }
    }
    else
    {
        forAll(tmpp, iface)
        {
            tmpValFrac[iface] = 1.0;
        }
    }
    fvPatchField<scalar>::operator=(tmpp);

    this->valueFraction() = tmpValFrac;
    this->refGrad() = poroHydraulic_->pField().dimensions()==dimLength ? scalarField(this->size(), 0.0) - n_ : scalarField(this->size(), 0.0);
    this->refValue() = tmpp;

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::timeDependentValveFvPatchScalarField::write(Ostream &os)
    const
{
    fvPatchScalarField::write(os);
    outletValue_.writeEntry("outletValue", os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        timeDependentValveFvPatchScalarField);
}

// ************************************************************************* //
