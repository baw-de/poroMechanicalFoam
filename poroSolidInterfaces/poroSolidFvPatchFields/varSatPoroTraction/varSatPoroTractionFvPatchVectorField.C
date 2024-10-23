/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "varSatPoroTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "poroFluidModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    effectiveStressModel& varSatPoroTractionFvPatchVectorField::effectiveStressModelRef()
    {
        if(!effectiveStressModelPtr_.valid())
        {
            makeEffectiveStressModel();
        }
        return effectiveStressModelPtr_.ref();
    }

    void varSatPoroTractionFvPatchVectorField::makeEffectiveStressModel()
    {
        dictionary dict;
        effectiveStressModelPtr_.reset(effectiveStressModel::New(dict,effectiveStressModelName_,this->patch().boundaryMesh().mesh()));
    }


    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    varSatPoroTractionFvPatchVectorField::
        varSatPoroTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF)
        : fixedGradientFvPatchVectorField(p, iF),
          totalTraction_(),
          traction_(p.size(), vector::zero),
          pressure_(p.size(), 0.0),
          tractionSeries_(),
          pressureSeries_(),
          effectiveStressModelName_(),
          effectiveStressModelPtr_(),
          secondOrder_(false),
          limitCoeff_(1.0),
          relaxFac_(1.0),
          buoyancyImplicit_(false)
    {
        fvPatchVectorField::operator=(patchInternalField());
        gradient() = vector::zero;
    }

    varSatPoroTractionFvPatchVectorField::
        varSatPoroTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const dictionary &dict)
        : fixedGradientFvPatchVectorField(p, iF),
          totalTraction_(dict.lookupOrDefault<Switch>("total", true)),
          traction_(p.size(), vector::zero),
          pressure_(p.size(), 0.0),
          tractionSeries_(PatchFunction1<vector>::New(p.patch(), "traction", dict)),
          pressureSeries_(PatchFunction1<scalar>::New(p.patch(), "pressure", dict)),
          effectiveStressModelName_(dict.lookupOrDefault<word>("effectiveStressModel","terzaghi")),
          effectiveStressModelPtr_(),
          secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
          limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0)),
          relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
          buoyancyImplicit_(false)
    {
        Info << "Creating " << type() << " boundary condition at patch " << patch().name() <<  nl;
        totalTraction_
        ? Info << tab << "Traction is total" << endl
        : Info << tab << "Traction is effective" << endl;

        effectiveStressModelPtr_.reset(effectiveStressModel::New(dict,effectiveStressModelName_,this->patch().boundaryMesh().mesh()));

        if (!totalTraction_)
        {
            buoyancyImplicit_ = dict.lookupOrDefault<Switch>("buoyancyIncluded", false);
            if(buoyancyImplicit_)
            {
                WarningInFunction() << "For variably saturated calculations, please dont include buoyancy in total traction"
                                    << " boundary conditions!" << endl;
            }
        }

        if (dict.found("gradient"))
        {
            gradient() = vectorField("gradient", dict, p.size());
        }
        else
        {
            gradient() = vector::zero;
        }

        if (dict.found("value"))
        {
            Field<vector>::operator=(vectorField("value", dict, p.size()));
        }
        else
        {
            fvPatchVectorField::operator=(patchInternalField());
        }

        if (secondOrder_)
        {
            Info << "    second order correction" << endl;
        }

        if (limitCoeff_)
        {
            Info << "    limiter coefficient: " << limitCoeff_ << endl;
        }

        if (relaxFac_ < 1.0)
        {
            Info << "    relaxation factor: " << relaxFac_ << endl;
        }
    }

    varSatPoroTractionFvPatchVectorField::
        varSatPoroTractionFvPatchVectorField(
            const varSatPoroTractionFvPatchVectorField &stpvf,
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const fvPatchFieldMapper &mapper)
        : fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
          totalTraction_(stpvf.totalTraction_),
          traction_(stpvf.traction_, mapper),
          pressure_(stpvf.pressure_, mapper),
          tractionSeries_(stpvf.tractionSeries_.clone(this->patch().patch())),
          pressureSeries_(stpvf.pressureSeries_.clone(this->patch().patch())),
          effectiveStressModelName_(stpvf.effectiveStressModelName_),
          effectiveStressModelPtr_(stpvf.effectiveStressModelPtr_),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_),
          buoyancyImplicit_(stpvf.buoyancyImplicit_)
    {
    }

    varSatPoroTractionFvPatchVectorField::
        varSatPoroTractionFvPatchVectorField(
            const varSatPoroTractionFvPatchVectorField &stpvf)
        : fixedGradientFvPatchVectorField(stpvf),
          totalTraction_(stpvf.totalTraction_),
          traction_(stpvf.traction_),
          pressure_(stpvf.pressure_),
          tractionSeries_(stpvf.tractionSeries_.clone(this->patch().patch())),
          pressureSeries_(stpvf.pressureSeries_.clone(this->patch().patch())),
          effectiveStressModelName_(stpvf.effectiveStressModelName_),
          effectiveStressModelPtr_(stpvf.effectiveStressModelPtr_),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_),
          buoyancyImplicit_(stpvf.buoyancyImplicit_)
    {
    }

    varSatPoroTractionFvPatchVectorField::
        varSatPoroTractionFvPatchVectorField(
            const varSatPoroTractionFvPatchVectorField &stpvf,
            const DimensionedField<vector, volMesh> &iF)
        : fixedGradientFvPatchVectorField(stpvf, iF),
          totalTraction_(stpvf.totalTraction_),
          traction_(stpvf.traction_),
          pressure_(stpvf.pressure_),
          tractionSeries_(stpvf.tractionSeries_.clone(this->patch().patch())),
          pressureSeries_(stpvf.pressureSeries_.clone(this->patch().patch())),
          effectiveStressModelName_(stpvf.effectiveStressModelName_),
          effectiveStressModelPtr_(stpvf.effectiveStressModelPtr_),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_),
          buoyancyImplicit_(stpvf.buoyancyImplicit_)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void varSatPoroTractionFvPatchVectorField::autoMap(
        const fvPatchFieldMapper &m)
    {
        fixedGradientFvPatchVectorField::autoMap(m);

        traction_.autoMap(m);
        pressure_.autoMap(m);
    }

    // Reverse-map the given fvPatchField onto this fvPatchField
    void varSatPoroTractionFvPatchVectorField::rmap(
        const fvPatchVectorField &ptf,
        const labelList &addr)
    {
        fixedGradientFvPatchVectorField::rmap(ptf, addr);

        const varSatPoroTractionFvPatchVectorField &dmptf =
            refCast<const varSatPoroTractionFvPatchVectorField>(ptf);

        traction_.rmap(dmptf.traction_, addr);
        pressure_.rmap(dmptf.pressure_, addr);
    }

    // Update the coefficients associated with the patch field
    void varSatPoroTractionFvPatchVectorField::updateCoeffs()
    {

        if (updated()||this->db().time().timeIndex()<this->db().time().startTimeIndex()+1)
        {
            return;
        }

        traction_ = tractionSeries_->value(this->db().time().timeOutputValue());
        pressure_ = pressureSeries_->value(this->db().time().timeOutputValue());


        scalarField totalPressure(pressure_);

     if(!totalTraction_)
     {
        //const fvPatchField<scalar> &p =
        const fvPatchField<scalar> p =
                    this->db().lookupObject<volScalarField>("p").boundaryField()[patch().index()];
        const fvPatchField<scalar> &S =
                    this->db().lookupObject<volScalarField>("S").boundaryField()[patch().index()];
        //const fvPatchField<scalar> &n =
        //                this->patch().patchField<volScalarField, scalar>(this->db().lookupObject<volScalarField>("n"));

        scalarField porosity(patch().size(),1.0);

        scalarField chi = effectiveStressModelRef().chi(porosity,S,p);
                totalPressure = pressure_ + (chi) * p; //- chi_ * p_;

        Info << "The prescribed (max:min) effective pressure of: " << max(pressure_) << ":" << min(pressure_) << " corresponds to " << nl
             << max(totalPressure) << ":" << min(totalPressure) << " in total pressure" << endl;
        //Info << min(mag(chi))
     }
        // Lookup the solidModel object
        const solidModel &solMod = lookupSolidModel(patch().boundaryMesh().mesh());

        // Set surface-normal gradient on the patch corresponding to the desired
        // traction
        gradient() =
            relaxFac_ * solMod.tractionBoundarySnGrad(
                            traction_, totalPressure, patch()) +
            (1.0 - relaxFac_) * gradient();

        fixedGradientFvPatchVectorField::updateCoeffs();
    }

    void varSatPoroTractionFvPatchVectorField::evaluate(
        const Pstream::commsTypes commsType)
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }

        // Lookup the gradient field
        const fvPatchField<tensor> &gradField =
            patch().lookupPatchField<volTensorField, tensor>(
                "grad(" + internalField().name() + ")"
            );

        // Face unit normals
        const vectorField n = patch().nf();

        // Delta vectors
        const vectorField delta = patch().delta();

        // Non-orthogonal correction vectors
        const vectorField k = ((I - sqr(n)) & delta);

        if (secondOrder_)
        {
            const vectorField dUP = (k & gradField.patchInternalField());
            const vectorField nGradUP = (n & gradField.patchInternalField());

            Field<vector>::operator=(
                patchInternalField() + dUP + 0.5 * (gradient() + nGradUP) / patch().deltaCoeffs());
        }
        else
        {

            Field<vector>::operator=(
                patchInternalField() + (k & gradField.patchInternalField()) + gradient() / patch().deltaCoeffs());
        }

        fvPatchField<vector>::evaluate();
    }

    void varSatPoroTractionFvPatchVectorField::write(Ostream &os) const
    {
        // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
        // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
        // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
        // fixedGradientFvPatchVectorField::write(os);
        fvPatchVectorField::write(os);



        tractionSeries_->writeData(os);

        pressureSeries_->writeData(os);

        os.writeKeyword("total")
            << totalTraction_ << token::END_STATEMENT << nl;
        os.writeKeyword("secondOrder")
            << secondOrder_ << token::END_STATEMENT << nl;
        os.writeKeyword("limitCoeff")
            << limitCoeff_ << token::END_STATEMENT << nl;
        os.writeKeyword("relaxationFactor")
            << relaxFac_ << token::END_STATEMENT << nl;
        if(totalTraction_)
        {
        os.writeKeyword("buoyancyImplicit")
            << buoyancyImplicit_ << token::END_STATEMENT << nl;
        }

        writeEntry("value", os);

    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField(fvPatchVectorField, varSatPoroTractionFvPatchVectorField);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
