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

\*---------------------------------------------------------------------------*/

#include "fixedPoroFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "fvPatch.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedPoroFluxFvPatchScalarField::fixedPoroFluxFvPatchScalarField(
    const fvPatch &p,
    const DimensionedField<scalar, volMesh> &iF)
    : fixedGradientFvPatchScalarField(p, iF),
      flux_(p.size(), pTraits<scalar>::zero),
      fluxSeries_(),
      isHead_(false),
      HMCoupled(dynamicFvMesh::defaultRegion)
{
}

Foam::fixedPoroFluxFvPatchScalarField::fixedPoroFluxFvPatchScalarField(
    const fvPatch &p,
    const DimensionedField<scalar, volMesh> &iF,
    const dictionary &dict)
    : fixedGradientFvPatchScalarField(p, iF),
      flux_(p.size(), 0.0),
      fluxSeries_(PatchFunction1<scalar>::New(p.patch(), "flux", dict)),
      isHead_(false),
      HMCoupled(dynamicFvMesh::defaultRegion)
{
    if(iF.dimensions()==dimLength)
    {isHead_ = true;}
    if(this->db().name()!=dynamicFvMesh::defaultRegion)
    {
        HMCoupled = "poroFluid";
    }
    flux_=fluxSeries_->value(this->db().time().timeOutputValue());
    gradient() = -(0) * flux_;
    fixedGradientFvPatchScalarField::updateCoeffs();
    fixedGradientFvPatchScalarField::evaluate();
}

Foam::fixedPoroFluxFvPatchScalarField::fixedPoroFluxFvPatchScalarField(
    const fixedPoroFluxFvPatchScalarField &ptf,
    const fvPatch &p,
    const DimensionedField<scalar, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
      flux_(ptf.flux_),
      fluxSeries_(ptf.fluxSeries_.clone(this->patch().patch())),
      isHead_(ptf.isHead_),
      HMCoupled(ptf.HMCoupled)
{
}

Foam::fixedPoroFluxFvPatchScalarField::fixedPoroFluxFvPatchScalarField(
    const fixedPoroFluxFvPatchScalarField &ptf)
    : fixedGradientFvPatchScalarField(ptf),
      flux_(ptf.flux_),
      fluxSeries_(ptf.fluxSeries_.clone(this->patch().patch())),
      isHead_(ptf.isHead_),
      HMCoupled(ptf.HMCoupled)
{
}

Foam::fixedPoroFluxFvPatchScalarField::fixedPoroFluxFvPatchScalarField(
    const fixedPoroFluxFvPatchScalarField &ptf,
    const DimensionedField<scalar, volMesh> &iF)
    : fixedGradientFvPatchScalarField(ptf, iF),
      flux_(ptf.flux_),
      fluxSeries_(ptf.fluxSeries_.clone(this->patch().patch())),
      isHead_(ptf.isHead_),
      HMCoupled(ptf.HMCoupled)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedPoroFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

        flux_ = fluxSeries_->value(this->db().time().timeOutputValue());
        // Info << "h at patch " << patch().name() << " = "
        //     << headSeries_(this->db().time().timeOutputValue())
        //      << endl;

if(isHead_)
{
    const scalarField& k_eff_ = 
        this->patch().patchField<surfaceScalarField, scalar>(
            this->db().time().subRegistry(HMCoupled).lookupObject<surfaceScalarField>("kEfff")
            );
    const UniformDimensionedField<vector> &gamma = this->db().time().subRegistry(HMCoupled).lookupObject<UniformDimensionedField<vector>>("gamma_water");
    scalarField n_ = patch().nf() & vector(gamma.value()).normalise();
    gradient() = (flux_) / max(k_eff_,SMALL) + n_;
}
else
{
    const scalarField& k_eff_ = 
        this->patch().patchField<surfaceScalarField, scalar>(this->db().time().subRegistry(HMCoupled).lookupObject<surfaceScalarField>("kEffbyGammaf"));
    gradient() = (flux_ / k_eff_);
}

    //Info << "Flux on Patch " << patch().name() << ": " << flux_ << " gradZ: " << n_ << " Gradient: " << (flux_/(k_*kr_))-n_ << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void Foam::fixedPoroFluxFvPatchScalarField::write(Ostream &os) const
{
    fvPatchScalarField::write(os);
    fluxSeries_().writeData(os);
    writeEntry("value", os);
    //gradient().writeEntry("gradient", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        fixedPoroFluxFvPatchScalarField);
}

// ************************************************************************* //
