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

#include "limitedHeadInfiltrationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::limitedHeadInfiltrationFvPatchScalarField::
    limitedHeadInfiltrationFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(p, iF),
      pName_("p_rgh"),
      kEffname_("kEffbyGammaf"),
      updatedTime_(-1.0),
      pMax_(p.size(), 0.0),
      fluxSeries_(),
      isHead_(false),
      HMCoupled(dynamicFvMesh::defaultRegion),
      explicit_(false)
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}

Foam::limitedHeadInfiltrationFvPatchScalarField::
    limitedHeadInfiltrationFvPatchScalarField(
        const limitedHeadInfiltrationFvPatchScalarField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : mixedFvPatchScalarField(ptf, p, iF, mapper),
      pName_(ptf.pName_),
      kEffname_(ptf.kEffname_),
      updatedTime_(ptf.updatedTime_),
      pMax_(ptf.pMax_, mapper),
      fluxSeries_(ptf.fluxSeries_.clone(this->patch().patch())),
      isHead_(ptf.isHead_),
      HMCoupled(ptf.HMCoupled),
      explicit_(ptf.explicit_)
{
}

Foam::limitedHeadInfiltrationFvPatchScalarField::
    limitedHeadInfiltrationFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
    : mixedFvPatchScalarField(p, iF),
      pName_(dict.lookupOrDefault<word>("pressureField", "p_rgh")),
      kEffname_(dict.lookupOrDefault<word>("k_eff", "kEffbyGammaf")),
      updatedTime_(-1.0),
      pMax_("pMax", dict, p.size()),
      fluxSeries_(PatchFunction1<scalar>::New(p.patch(), "flux", dict)),
      isHead_(false),
      HMCoupled(dynamicFvMesh::defaultRegion),
      explicit_(dict.lookupOrDefault<Switch>("explicit",false))
{
    if(iF.dimensions()==dimLength)
    {
        isHead_ = true;
        pName_ = "pHead";
        kEffname_ = "kEfff";
    }
    if(!(this->db().name()==dynamicFvMesh::defaultRegion))
    {
        HMCoupled = "poroFluid";
    }
    this->refValue() = pTraits<scalar>::zero;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=(
            scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchField<scalar>::operator=(updatedTime_);
    }

    this->refGrad() = scalarField(this->size(), 0.0);
    this->valueFraction() = 0.0;
}

Foam::limitedHeadInfiltrationFvPatchScalarField::
    limitedHeadInfiltrationFvPatchScalarField(
        const limitedHeadInfiltrationFvPatchScalarField &tppsf)
    : mixedFvPatchScalarField(tppsf),
      pName_(tppsf.pName_),
      kEffname_(tppsf.kEffname_),
      updatedTime_(tppsf.updatedTime_),
      pMax_(tppsf.pMax_),
      fluxSeries_(tppsf.fluxSeries_.clone(this->patch().patch())),
      isHead_(tppsf.isHead_),
      HMCoupled(tppsf.HMCoupled),
      explicit_(tppsf.explicit_)
{
}

Foam::limitedHeadInfiltrationFvPatchScalarField::
    limitedHeadInfiltrationFvPatchScalarField(
        const limitedHeadInfiltrationFvPatchScalarField &tppsf,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(tppsf, iF),
      pName_(tppsf.pName_),
      kEffname_(tppsf.kEffname_),
      updatedTime_(tppsf.updatedTime_),
      pMax_(tppsf.pMax_),
      fluxSeries_(tppsf.fluxSeries_.clone(this->patch().patch())),
      isHead_(tppsf.isHead_),
      HMCoupled(tppsf.HMCoupled),
      explicit_(tppsf.explicit_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::limitedHeadInfiltrationFvPatchScalarField::autoMap(
    const fvPatchFieldMapper &m)
{
    mixedFvPatchScalarField::autoMap(m);
    pMax_.autoMap(m);
}

void Foam::limitedHeadInfiltrationFvPatchScalarField::rmap(
    const fvPatchScalarField &ptf,
    const labelList &addr)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const limitedHeadInfiltrationFvPatchScalarField &tiptf =
        refCast<const limitedHeadInfiltrationFvPatchScalarField>(ptf);

    pMax_.rmap(tiptf.pMax_, addr);
}

void Foam::limitedHeadInfiltrationFvPatchScalarField::evaluate(
    const Pstream::commsTypes commsType)
{
    mixedFvPatchScalarField::evaluate();
}

void Foam::limitedHeadInfiltrationFvPatchScalarField::updateCoeffs()
{
   
    if (updated())
    {
        return;
    }

    bool allowUpdate = false;

    if(!explicit_ || updatedTime_ < this->db().time().timeOutputValue())
    {
        allowUpdate = true;
        updatedTime_ = this->db().time().timeOutputValue();
        //Info << "Updating boundary " << patch().name() << " with type " << typeName << explicit_?" implicitly":" explicitly" << endl;
    }

    if(allowUpdate)
    { 
	    scalarField flux_ = fluxSeries_->value(this->db().time().timeOutputValue());

	    const fvsPatchField<scalar> &kEff_ =
		this->patch().patchField<surfaceScalarField, scalar>(this->db().lookupObject<surfaceScalarField>(kEffname_));

	    
	    scalarField wantedGrad(this->size(), 0.0);
	    scalarField wantedVal(this->size(), 0.0);

	    if(isHead_)
	    {
		    const UniformDimensionedField<vector> &gamma = this->db().lookupObject<UniformDimensionedField<vector>>("gamma_water");
		    scalarField n_ = patch().nf() & vector(gamma.value()).normalise();
		    wantedGrad = ((flux_ / (kEff_)) + n_);
		    wantedVal = pMax_;
	    }
	    else
	    {
		    const fvPatchField<scalar> &pHyd =
			this->patch().patchField<volScalarField, scalar>(this->db().lookupObject<volScalarField>("p_Hyd"));
		    wantedGrad = (flux_ / (kEff_));
		    wantedVal = pMax_ - pHyd;
	    }
	    
	    scalarField tmpValFrac(this->size(), 0.0);
	    	    
	    if(isHead_)
	    {	
	    	const fvPatchField<scalar> &pField =
		    this->patch().patchField<volScalarField, scalar>(this->db().lookupObject<volScalarField>("pHead"));
	    	tmpValFrac = pos(pField - pMax_);
	    }
	    else
	    {
	        const fvPatchField<scalar> &pTotal =
			this->patch().patchField<volScalarField, scalar>(this->db().lookupObject<volScalarField>("p"));
	    	tmpValFrac = pos(pTotal - pMax_);
	    }
	    	    
	    const scalarField &currentGrad = snGrad();
	    forAll(wantedGrad, iFace)
		{
			if (currentGrad[iFace] > wantedGrad[iFace])
			{
			    tmpValFrac[iFace] = 0;
			}
        }   
		    
		this->valueFraction() = tmpValFrac;
	

        this->refGrad() = wantedGrad;
        this->refValue() = wantedVal;

        mixedFvPatchScalarField::updateCoeffs();
    }
}

void Foam::limitedHeadInfiltrationFvPatchScalarField::write(Ostream &os)
    const
{
    fvPatchScalarField::write(os);
    //os.writeKeyword("flux") << updatedTime_ << token::END_STATEMENT << nl;
    fluxSeries_().writeData(os);
    pMax_.writeEntry("pMax", os);
    os.writeEntry("pressureField", pName_);
    os.writeEntry("k_eff", kEffname_);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        limitedHeadInfiltrationFvPatchScalarField);
}

// ************************************************************************* //
