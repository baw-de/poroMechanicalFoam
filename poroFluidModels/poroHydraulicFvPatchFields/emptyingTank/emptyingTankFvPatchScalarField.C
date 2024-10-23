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

#include "emptyingTankFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emptyingTankFvPatchScalarField::
    emptyingTankFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(p, iF),
      pName_("pHead"),
      fluxname_("phi"),
      h0_(0.0),
      dh_(0.0),
      time_(patch().boundaryMesh().mesh().time()),
      crossSectionSeries_(),
      HMCoupled(dynamicFvMesh::defaultRegion)
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}

Foam::emptyingTankFvPatchScalarField::
    emptyingTankFvPatchScalarField(
        const emptyingTankFvPatchScalarField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : mixedFvPatchScalarField(ptf, p, iF, mapper),
      pName_(ptf.pName_),
      fluxname_(ptf.fluxname_),
      h0_(ptf.h0_),
      dh_(ptf.dh_),
      time_(ptf.time_),
      crossSectionSeries_(ptf.crossSectionSeries_),
      HMCoupled(ptf.HMCoupled)
{}

Foam::emptyingTankFvPatchScalarField::
    emptyingTankFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
    : mixedFvPatchScalarField(p, iF),
      pName_(dict.lookupOrDefault<word>("pField", "p_rgh")),
      fluxname_(dict.lookupOrDefault<word>("flux", "phi")),
      h0_(readScalar(dict.lookup("h0"))),
      dh_(0.0),
      time_(patch().db().time()),
      crossSectionSeries_(),
      HMCoupled(dynamicFvMesh::defaultRegion)
{
    if(!(this->db().name()==dynamicFvMesh::defaultRegion))
    {
        HMCoupled = "poroFluid";
    }

    this->refValue() = pTraits<scalar>::zero;

    if (dict.found("crossSectionSeries"))
    {
        Info << " tank cross section is non-uniform" << endl;
        crossSectionSeries_ =
            interpolationTable<scalar>(dict.subDict("crossSectionSeries"));
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=(
            scalarField("value", dict, p.size()));
    }

    this->refGrad() = scalarField(this->size(), 0.0);
    this->valueFraction() = 0.0;
}

Foam::emptyingTankFvPatchScalarField::
    emptyingTankFvPatchScalarField(
        const emptyingTankFvPatchScalarField &tppsf)
    : mixedFvPatchScalarField(tppsf),
      pName_(tppsf.pName_),
      fluxname_(tppsf.fluxname_),
      h0_(tppsf.h0_),
      dh_(tppsf.dh_),
      time_(tppsf.time_),
      crossSectionSeries_(tppsf.crossSectionSeries_),
      HMCoupled(tppsf.HMCoupled)
{
}

Foam::emptyingTankFvPatchScalarField::
    emptyingTankFvPatchScalarField(
        const emptyingTankFvPatchScalarField &tppsf,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(tppsf, iF),
      pName_(tppsf.pName_),
      fluxname_(tppsf.fluxname_),
      h0_(tppsf.h0_),
      dh_(tppsf.dh_),
      time_(tppsf.time_),
      crossSectionSeries_(tppsf.crossSectionSeries_),
      HMCoupled(tppsf.HMCoupled)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::emptyingTankFvPatchScalarField::autoMap(
    const fvPatchFieldMapper &m)
{
    mixedFvPatchScalarField::autoMap(m);
}

void Foam::emptyingTankFvPatchScalarField::rmap(
    const fvPatchScalarField &ptf,
    const labelList &addr)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}

void Foam::emptyingTankFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const poroFluidModel& pFM = db().time().lookupObject<poroFluidModel>(HMCoupled);

    dimensionedScalar currTime = this->db().time();
    if (time_.value() != currTime.value())
    {
        h0_ = h0_ + dh_;
        time_ = currTime;
    }

    scalar crossSection_ = 1.0;
    if (crossSectionSeries_.size())
    {
        crossSection_ = crossSectionSeries_(this->db().time().timeOutputValue());
    }

    const fvsPatchField<scalar> &flux =
        this->patch().patchField<surfaceScalarField, scalar>(pFM.phi());
    
        scalarField n_ = patch().nf() & vector(pFM.poroHydraulic().gamma().value()).normalise();
        scalarField z_ = patch().Cf() & vector(pFM.poroHydraulic().gamma().value()).normalise();

    scalarField tmpValFrac(this->size(), 1.0);
    dh_ = gSum(flux) / (crossSection_ * gSum(patch().magSf())) * patch().boundaryMesh().mesh().time().deltaT().value();

    forAll(tmpValFrac, iFace)
    {
        if (h0_ + 0.5 * dh_ < z_[iFace])
        {
            tmpValFrac[iFace] = 0;
        }
        else
        {
            tmpValFrac[iFace] = 1;
        }
    }

    this->valueFraction() = tmpValFrac;
    if(pFM.pField().dimensions() == dimLength){
    this->refGrad() = (-n_);
    this->refValue() = h0_ + 0.5 * dh_ - z_;
    }
    else{
    this->refGrad() = scalarField(this->size(),0.0);
    this->refValue() = (h0_ + 0.5 * dh_ - pFM.poroHydraulic().href().value())*pFM.poroHydraulic().gamma().value();
    }

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::emptyingTankFvPatchScalarField::write(Ostream &os)
    const
{
    fvPatchScalarField::write(os);
    os.writeEntry("h0",h0_);
    os.writeEntry("pField",pName_);
    os.writeEntry("flux",fluxname_);
    if(crossSectionSeries_.valid())
    {
        crossSectionSeries_().writeData(os);
    }
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        emptyingTankFvPatchScalarField);
}

// ************************************************************************* //
