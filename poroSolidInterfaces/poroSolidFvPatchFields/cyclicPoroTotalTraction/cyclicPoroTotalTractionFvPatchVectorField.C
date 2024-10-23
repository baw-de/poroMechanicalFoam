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

#include "cyclicPoroTotalTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    cyclicPoroTotalTractionFvPatchVectorField::
        cyclicPoroTotalTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF)
        : fixedGradientFvPatchVectorField(p, iF),
          tractionAmp_(p.size(), vector::zero),
          pressureAmp_(p.size(), 0.0),
          tractionBase_(p.size(), vector::zero),
          pressureBase_(p.size(), 0.0),
          T_(0.0),
          f_(0.0),
          shift_(0.0),
          secondOrder_(false),
          limitCoeff_(1.0),
          relaxFac_(1.0)
    {
        fvPatchVectorField::operator=(patchInternalField());
        gradient() = vector::zero;
    }

    cyclicPoroTotalTractionFvPatchVectorField::
        cyclicPoroTotalTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const dictionary &dict)
        : fixedGradientFvPatchVectorField(p, iF),
          tractionAmp_(p.size(), vector::zero),
          pressureAmp_(p.size(), 0.0),
          tractionBase_(p.size(), vector::zero),
          pressureBase_(p.size(), 0.0),
          T_(1.0),
          f_(1.0),
          shift_(dict.lookupOrDefault<scalar>("phaseShift", 0.0)),
          secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
          limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0)),
          relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))
    {
        Info << "Creating " << type() << " boundary condition" << endl;

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

        tractionAmp_ = vectorField("tractionAmplitude", dict, p.size());
        tractionBase_ = vectorField("BaseTraction", dict, p.size());
        pressureAmp_ = scalarField("pressureAmplitude", dict, p.size());
        pressureBase_ = scalarField("BasePressure", dict, p.size());

        if (dict.found("period"))
        {
            T_ = scalar(readScalar(dict.lookup("period")));
        }
        else
        {
            f_ = scalar(readScalar(dict.lookup("frequency")));
            T_ = 1 / f_;
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

    cyclicPoroTotalTractionFvPatchVectorField::
        cyclicPoroTotalTractionFvPatchVectorField(
            const cyclicPoroTotalTractionFvPatchVectorField &stpvf,
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const fvPatchFieldMapper &mapper)
        : fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
          tractionAmp_(mapper(stpvf.tractionAmp_)),
          pressureAmp_(mapper(stpvf.pressureAmp_)),
          tractionBase_(mapper(stpvf.tractionBase_)),
          pressureBase_(mapper(stpvf.pressureBase_)),
#else
          tractionAmp_(stpvf.tractionAmp_, mapper),
          pressureAmp_(stpvf.pressureAmp_, mapper),
          tractionBase_(stpvf.tractionBase_, mapper),
          pressureBase_(stpvf.pressureBase_, mapper),
#endif
          T_(stpvf.T_),
          f_(stpvf.f_),
          shift_(stpvf.shift_),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_)
    {
    }

    cyclicPoroTotalTractionFvPatchVectorField::
        cyclicPoroTotalTractionFvPatchVectorField(
            const cyclicPoroTotalTractionFvPatchVectorField &stpvf)
        : fixedGradientFvPatchVectorField(stpvf),
          tractionAmp_(stpvf.tractionAmp_),
          pressureAmp_(stpvf.pressureAmp_),
          tractionBase_(stpvf.tractionBase_),
          pressureBase_(stpvf.pressureBase_),
          T_(stpvf.T_),
          f_(stpvf.f_),
          shift_(stpvf.shift_),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_)
    {
    }

    cyclicPoroTotalTractionFvPatchVectorField::
        cyclicPoroTotalTractionFvPatchVectorField(
            const cyclicPoroTotalTractionFvPatchVectorField &stpvf,
            const DimensionedField<vector, volMesh> &iF)
        : fixedGradientFvPatchVectorField(stpvf, iF),
          tractionAmp_(stpvf.tractionAmp_),
          pressureAmp_(stpvf.pressureAmp_),
          tractionBase_(stpvf.tractionBase_),
          pressureBase_(stpvf.pressureBase_),
          T_(stpvf.T_),
          f_(stpvf.f_),
          shift_(stpvf.shift_),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void cyclicPoroTotalTractionFvPatchVectorField::autoMap(
        const fvPatchFieldMapper &m)
    {
        fixedGradientFvPatchVectorField::autoMap(m);
#ifdef OPENFOAMFOUNDATION
        m(tractionAmp_, tractionAmp_);
        m(pressureAmp_, pressureAmp_);
        m(tractionBase_, tractionBase_);
        m(pressureBase_, pressureBase_);
#else
        tractionAmp_.autoMap(m);
        pressureAmp_.autoMap(m);
        tractionBase_.autoMap(m);
        pressureBase_.autoMap(m);
#endif
    }

    // Reverse-map the given fvPatchField onto this fvPatchField
    void cyclicPoroTotalTractionFvPatchVectorField::rmap(
        const fvPatchVectorField &ptf,
        const labelList &addr)
    {
        fixedGradientFvPatchVectorField::rmap(ptf, addr);

        const cyclicPoroTotalTractionFvPatchVectorField &dmptf =
            refCast<const cyclicPoroTotalTractionFvPatchVectorField>(ptf);

        tractionAmp_.rmap(dmptf.tractionAmp_, addr);
        pressureAmp_.rmap(dmptf.pressureAmp_, addr);
        tractionBase_.rmap(dmptf.tractionBase_, addr);
        pressureBase_.rmap(dmptf.pressureBase_, addr);
    }

    // Update the coefficients associated with the patch field
    void cyclicPoroTotalTractionFvPatchVectorField::updateCoeffs()
    {
        if (updated())
        {
            return;
        }

        const fvPatchField<scalar> &p_ =
            patch().lookupPatchField<volScalarField, scalar>("pHead");
        const fvPatchField<scalar> &chi_ =
            patch().lookupPatchField<volScalarField, scalar>("chi");
        const uniformDimensionedVectorField &gamma_w_ =
            db().lookupObject<uniformDimensionedVectorField>("gamma_w");

        const scalar t = db().time().value();

        scalar omega = 2 * acos(-1.0) / T_;

        scalarField effPressure = -(pressureBase_ + pressureAmp_ * sin(omega * t + shift_)) - chi_ * p_ * mag(gamma_w_.value());
        vectorField traction_ = tractionBase_ + tractionAmp_ * sin(omega * t + shift_);

        // Lookup the solidModel object
        const solidModel &solMod = lookupSolidModel(patch().boundaryMesh().mesh());

        // Set surface-normal gradient on the patch corresponding to the desired
        // traction
        gradient() =
            relaxFac_ * solMod.tractionBoundarySnGrad(
                            traction_, effPressure, patch()) +
            (1.0 - relaxFac_) * gradient();

        fixedGradientFvPatchVectorField::updateCoeffs();
    }

    void cyclicPoroTotalTractionFvPatchVectorField::evaluate(
        const Pstream::commsTypes commsType)
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }

        // Lookup the gradient field
        const fvPatchField<tensor> &gradField =
            patch().lookupPatchField<volTensorField, tensor>(
#ifdef OPENFOAMESIORFOUNDATION
                "grad(" + internalField().name() + ")"
#else
                "grad(" + dimensionedInternalField().name() + ")"
#endif
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

    void cyclicPoroTotalTractionFvPatchVectorField::write(Ostream &os) const
    {
        // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
        // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
        // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
        //fixedGradientFvPatchVectorField::write(os);
        fvPatchVectorField::write(os);

#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "tractionAmplitude", tractionAmp_);
        writeEntry(os, "BaseTraction", tractionBase_);
#else
        tractionAmp_.writeEntry("tractionAmplitude", os);
        tractionBase_.writeEntry("BaseTraction", os);
#endif

#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "pressureAmplitude", pressureAmp_);
        writeEntry(os, "BasePressure", pressureBase_);
#else
        pressureAmp_.writeEntry("pressureAmplitude", os);
        pressureBase_.writeEntry("BasePressure", os);
#endif
        os.writeKeyword("period")
            << T_ << token::END_STATEMENT << nl;
        os.writeKeyword("phaseShift")
            << shift_ << token::END_STATEMENT << nl;
        os.writeKeyword("secondOrder")
            << secondOrder_ << token::END_STATEMENT << nl;
        os.writeKeyword("limitCoeff")
            << limitCoeff_ << token::END_STATEMENT << nl;
        os.writeKeyword("relaxationFactor")
            << relaxFac_ << token::END_STATEMENT << nl;

#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "value", *this);
#else
        writeEntry("value", os);
#endif
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField(fvPatchVectorField, cyclicPoroTotalTractionFvPatchVectorField);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
