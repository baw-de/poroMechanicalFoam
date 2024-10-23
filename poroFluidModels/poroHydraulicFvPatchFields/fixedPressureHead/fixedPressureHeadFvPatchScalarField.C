/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fixedPressureHeadFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedPressureHeadFvPatchScalarField::fixedPressureHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(p.size(), 0.0),
    pressureSeries_(),
    isHead_(),
    HMCoupled(dynamicFvMesh::defaultRegion)
{}


Foam::fixedPressureHeadFvPatchScalarField::fixedPressureHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(p.size(), 0.0),
    pressureSeries_(PatchFunction1<scalar>::New(p.patch(), "pHead", dict)),
    isHead_(false),
    HMCoupled(dynamicFvMesh::defaultRegion)
{


    if(iF.dimensions()==dimLength)
    {isHead_ = true;}
    if(!(this->db().name()==dynamicFvMesh::defaultRegion))
    {
        HMCoupled = "poroFluid";
    }
    Info << "Creating porohydraulic boundary condition:  prescribed pressure head:" << endl;
    if (dict.found("pressureSeries"))
    {
        Info<< "  Pressure head at boundary " << patch().name() << " is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else {
      Info << "   Pressure head at boundary " << patch().name() << " is constant" << endl;
        p0_ = scalarField("pHead", dict, p.size());
    }

    //poroHydraulic_ = patch().db().lookupObject<poroHydraulicModel>("poroHydraulicModel");

    if (dict.found("value"))
      {
          fvPatchField<scalar>::operator=
          (
              scalarField("value", dict, p.size())
          );
      }
    else if (isHead_) {
          fvPatchField<scalar>::operator=
          (
              p0_
          );
      }
    else {
        const poroHydraulicModel* poroHydraulic_ = &this->db().time().lookupObject<poroHydraulicModel>("poroHydraulicModel");
        const fvPatchField<scalar>& z =
            patch().patchField<volScalarField, scalar>(poroHydraulic_->z());
        fvPatchField<scalar>::operator=
          (
              (p0_ + z - poroHydraulic_->href().value())*poroHydraulic_->magGamma().value()
          );
    }
}


Foam::fixedPressureHeadFvPatchScalarField::fixedPressureHeadFvPatchScalarField
(
    const fixedPressureHeadFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_, mapper),
    pressureSeries_(ptf.pressureSeries_.clone(this->patch().patch())),
    isHead_(ptf.isHead_),
    HMCoupled(ptf.HMCoupled)
{}


Foam::fixedPressureHeadFvPatchScalarField::fixedPressureHeadFvPatchScalarField
(
    const fixedPressureHeadFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    p0_(tppsf.p0_),
    pressureSeries_(tppsf.pressureSeries_.clone(this->patch().patch())),
    isHead_(tppsf.isHead_),
    HMCoupled(tppsf.HMCoupled)
{}


Foam::fixedPressureHeadFvPatchScalarField::fixedPressureHeadFvPatchScalarField
(
    const fixedPressureHeadFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    p0_(tppsf.p0_),
    pressureSeries_(tppsf.pressureSeries_.clone(this->patch().patch())),
    isHead_(tppsf.isHead_),
    HMCoupled(tppsf.HMCoupled)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedPressureHeadFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}


void Foam::fixedPressureHeadFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fixedPressureHeadFvPatchScalarField& tiptf =
        refCast<const fixedPressureHeadFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::fixedPressureHeadFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

        
    	p0_ = pressureSeries_->value(this->db().time().timeOutputValue());
        //Info << "h at patch " << patch().name() << " = "
        //    << pressureSeries_(this->db().time().timeOutputValue())
        //     << endl;

	if (isHead_)
	{
	    operator==
	    (
		p0_
	    );
	}
	else
	{
        const poroHydraulicModel *poroHydraulic_ = &this->db().time().lookupObject<poroHydraulicModel>("poroHydraulicModel");
	    const fvPatchField<scalar>& z =
		    patch().patchField<volScalarField, scalar>(poroHydraulic_->z());
	    operator==
	    (
		(p0_ + z - poroHydraulic_->href().value())*poroHydraulic_->magGamma().value()
	    );
	}

    fixedValueFvPatchScalarField::updateCoeffs();
}



void Foam::fixedPressureHeadFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    p0_.writeEntry("pHead", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedPressureHeadFvPatchScalarField
    );
}

// ************************************************************************* //
