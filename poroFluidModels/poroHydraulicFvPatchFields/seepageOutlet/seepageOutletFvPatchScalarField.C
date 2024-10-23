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

#include "seepageOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seepageOutletFvPatchScalarField::
    seepageOutletFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(p, iF),
      outletPressure_(0.0),
      h0_(p.size(), 0.0),
      hSeries_(),
      depuitApprox_(false),
      isHead_(false),
      HMCoupled(dynamicFvMesh::defaultRegion),
      explicit_(false),
      updatedTime_(-1.0)
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}

Foam::seepageOutletFvPatchScalarField::
    seepageOutletFvPatchScalarField(
        const seepageOutletFvPatchScalarField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : mixedFvPatchScalarField(ptf, p, iF, mapper),
      outletPressure_(ptf.outletPressure_),
      h0_(ptf.h0_, mapper),
      hSeries_(ptf.hSeries_),
      depuitApprox_(ptf.depuitApprox_),
      isHead_(ptf.isHead_),
      HMCoupled(ptf.HMCoupled),
      explicit_(ptf.explicit_),
      updatedTime_(ptf.updatedTime_)
{
}

Foam::seepageOutletFvPatchScalarField::
    seepageOutletFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
    : mixedFvPatchScalarField(p, iF),
      outletPressure_("outletPressure", dict, p.size()),
      h0_(p.size(), 0.0),
      hSeries_(PatchFunction1<scalar>::New(p.patch(), "h", dict)),
      depuitApprox_(dict.lookupOrDefault("depuitApproximation", false)),
      isHead_(false),
      HMCoupled(dynamicFvMesh::defaultRegion),
      explicit_(dict.lookupOrDefault("explicit", false)),
      updatedTime_(-1.0)
{
    if(iF.dimensions()==dimLength)
    {isHead_ = true;}
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
        fvPatchField<scalar>::operator=(outletPressure_);
    }

    this->refGrad() = scalarField(this->size(), 0.0);
    this->valueFraction() = 0.0;
}

Foam::seepageOutletFvPatchScalarField::
    seepageOutletFvPatchScalarField(
        const seepageOutletFvPatchScalarField &tppsf)
    : mixedFvPatchScalarField(tppsf),
      outletPressure_(tppsf.outletPressure_),
      h0_(tppsf.h0_),
      hSeries_(tppsf.hSeries_.clone(this->patch().patch())),
      depuitApprox_(tppsf.depuitApprox_),
      isHead_(tppsf.isHead_),
      HMCoupled(tppsf.HMCoupled),
      explicit_(tppsf.explicit_),
      updatedTime_(tppsf.updatedTime_)
{
}

Foam::seepageOutletFvPatchScalarField::
    seepageOutletFvPatchScalarField(
        const seepageOutletFvPatchScalarField &tppsf,
        const DimensionedField<scalar, volMesh> &iF)
    : mixedFvPatchScalarField(tppsf, iF),
      outletPressure_(tppsf.outletPressure_),
      h0_(tppsf.h0_),
      hSeries_(tppsf.hSeries_.clone(this->patch().patch())),
      depuitApprox_(tppsf.depuitApprox_),
      isHead_(tppsf.isHead_),
      HMCoupled(tppsf.HMCoupled),
      explicit_(tppsf.explicit_),
      updatedTime_(tppsf.updatedTime_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::seepageOutletFvPatchScalarField::autoMap(
    const fvPatchFieldMapper &m)
{
    mixedFvPatchScalarField::autoMap(m);
    h0_.autoMap(m);
}

void Foam::seepageOutletFvPatchScalarField::rmap(
    const fvPatchScalarField &ptf,
    const labelList &addr)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const seepageOutletFvPatchScalarField &tiptf =
        refCast<const seepageOutletFvPatchScalarField>(ptf);

    h0_.rmap(tiptf.h0_, addr);
}

void Foam::seepageOutletFvPatchScalarField::updateCoeffs()
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
        Info << "Updating boundary " << patch().name() <<  endl;
    }

    if(allowUpdate)
    {
        h0_ = hSeries_->value(this->db().time().timeOutputValue());

        const poroHydraulicModel *poroHydraulic_ = &this->db().time().lookupObject<poroHydraulicModel>("poroHydraulicModel");
        const fvPatchField<scalar> &z_ =
                    this->patch().patchField<volScalarField, scalar>(poroHydraulic_->z());
        const fvsPatchField<scalar> &phiP_ =
                    this->patch().patchField<surfaceScalarField, scalar>
                    (
                        this->db().time().subRegistry(HMCoupled).lookupObject<surfaceScalarField>
                        (
                            "phi"
                        )
                    );

        scalar href = poroHydraulic_->href().value();

        const scalarField& thisField(*this);

        scalarField tmpp(this->size(), 0.0);
        scalarField tmpValFrac(this->size(), 0.0);
        if(isHead_)
        {
            // Set boundary pressure to 0 to calculate the ficticous gradient
            //this->refValue() = outletPressure_;
            // dz/dn
            //const scalarField n_ = patch().nf() & vector(poroHydraulic_->gamma().value()).normalise();
            //const scalarField& pif(this->patchInternalField());
            // get the ficticious gradient (dh/dn = dp/dn + dz/dn)
            //tmp<scalarField> toutFlowMarker(pos((this->refValue() - pif)*this->patch().deltaCoeffs()+n_));
            //const scalarField& outFlowMarker = toutFlowMarker();
            forAll(tmpp, iface)
            {
                // Is cellface under the external watertable?
                if (z_[iface] <= h0_[iface]+href)
                {
                    // fixed Value
                    tmpValFrac[iface] = 1.0;
                    // boundary value = hydrostatic pressure
                    tmpp[iface] = h0_[iface] - z_[iface];
                }
                // cellface is above the external watertable
                else
                {
                    // is the ficticous gradient outflowing?
                    if (thisField[iface] >= 0 && phiP_[iface] >= 0 && !depuitApprox_)
                    {
                        // fixed Value
                        tmpValFrac[iface] = 1.0;
                        // boundary value = outlet pressure (typically atmospheric)
                        tmpp[iface] = outletPressure_[iface];
                    }
                    // ficticious gradient inflowing
                    else
                    {
                        // fixed gradient
                        tmpValFrac[iface] = 0.0;
                        // boundary field value = internal value for cosmetics (zero normal gradient)
                        tmpp[iface] = patchInternalField()()[iface];
                    }
                }
            }
        }
        else
        {
            const fvPatchField<scalar> &pHydP_ =
                    this->patch().patchField<volScalarField, scalar>(poroHydraulic_->p_Hyd());

            const scalarField pP = thisField + pHydP_;
            // Set boundary pressure to 0 to calculate the ficticous gradient (p_rgh = 0 - pHyd)
            //this->refValue() = outletPressure_ - pHydP_;
            //const scalarField& pif(this->patchInternalField());
            // get the ficticious gradient
            //tmp<scalarField> toutFlowMarker(pos((this->refValue() - pif)*this->patch().deltaCoeffs()));
            //const scalarField& outFlowMarker = toutFlowMarker();

            //Info << outFlowMarker << endl;

            forAll(tmpp, iface)
            {
                // Is cellface under the external watertable?
                if (z_[iface] <= h0_[iface]+href)
                {
                    //Info << "bin drunter" << endl;
                    // fixed Value
                    tmpValFrac[iface] = 1.0;
                    // boundary value = watertable (i.e. potential) * gamma_w
                    tmpp[iface] =  h0_[iface]*poroHydraulic_->magGamma().value();
                }
                else
                {
                    // boundary value = outlet pressure (typically atmospheric) - hydrostratic pressure
                    tmpp[iface] = outletPressure_[iface] - pHydP_[iface];
                    if(this->valueFraction()[iface] == 0)
                    {
                        if (pP[iface]>outletPressure_[iface])
                        {
                            // fixed Value
                            tmpValFrac[iface] = 1.0;
                        }
                        else
                        {
                            // fixed gradient
                            tmpValFrac[iface] = 0.0;
                        }
                    }
                    else
                    {
                        if(phiP_[iface]>0)
                        {
                            // fixed gradient
                            tmpValFrac[iface] = 0.0;
                        }
                        else
                        {
                            // fixed Value
                            tmpValFrac[iface] = 1.0;
                        }
                    }
                }
            }
        }


        fvPatchField<scalar>::operator=(tmpp);

        this->valueFraction() = tmpValFrac;
        if(isHead_)
        {
            scalarField n_ = patch().nf() & vector(poroHydraulic_->gamma().value()).normalise();
            this->refGrad() = scalarField(this->size(), 0.0) - n_;
        }
        else
        {
            this->refGrad() = scalarField(this->size(), 0.0);
        }
        this->refValue() = tmpp;

        mixedFvPatchScalarField::updateCoeffs();
    }
}

void Foam::seepageOutletFvPatchScalarField::write(Ostream &os)
    const
{
    fvPatchScalarField::write(os);
    hSeries_().writeData(os);
    outletPressure_.writeEntry("outletPressure", os);

    os.writeEntry("depuitApproximation",depuitApprox_);
    os.writeEntry("explicit",explicit_);
    os.writeEntry("z",zName_);

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        seepageOutletFvPatchScalarField);
}

// ************************************************************************* //
