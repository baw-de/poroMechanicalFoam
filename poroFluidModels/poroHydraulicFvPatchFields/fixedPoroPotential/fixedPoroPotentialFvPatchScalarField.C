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

#include "fixedPoroPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedPoroPotentialFvPatchScalarField::fixedPoroPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    h0_(p.size(), 0.0),
    headSeries_(),
    isHead_(false),
    HMCoupled(dynamicFvMesh::defaultRegion)
{}


Foam::fixedPoroPotentialFvPatchScalarField::fixedPoroPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    h0_(p.size(), 0.0),
    headSeries_(PatchFunction1<scalar>::New(p.patch(), "h", dict)),
    isHead_(false),
    HMCoupled("region0")
{
    if(iF.dimensions()==dimLength)
    {isHead_ = true;}
    if(this->db().name()!=dynamicFvMesh::defaultRegion)
    {
        HMCoupled = "poroFluid";
    }
    Info << "Creating porohydraulic boundary condition:  prescribed total head:" << endl;

    //poroHydraulic_ = patch().db().lookupObject<poroHydraulicModel>("poroHydraulicModel");

    if (dict.found("value"))
      {
          fvPatchField<scalar>::operator=
          (
              scalarField("value", dict, p.size())
          );
      }
}


Foam::fixedPoroPotentialFvPatchScalarField::fixedPoroPotentialFvPatchScalarField
(
    const fixedPoroPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    h0_(ptf.h0_, mapper),
    headSeries_(ptf.headSeries_.clone(this->patch().patch())),
    isHead_(ptf.isHead_),
    HMCoupled(ptf.HMCoupled)
{}


Foam::fixedPoroPotentialFvPatchScalarField::fixedPoroPotentialFvPatchScalarField
(
    const fixedPoroPotentialFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    h0_(tppsf.h0_),
    headSeries_(tppsf.headSeries_.clone(this->patch().patch())),
    isHead_(tppsf.isHead_),
    HMCoupled(tppsf.HMCoupled)
{}


Foam::fixedPoroPotentialFvPatchScalarField::fixedPoroPotentialFvPatchScalarField
(
    const fixedPoroPotentialFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    h0_(tppsf.h0_),
    headSeries_(tppsf.headSeries_.clone(this->patch().patch())),
    isHead_(tppsf.isHead_),
    HMCoupled(tppsf.HMCoupled)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedPoroPotentialFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    h0_.autoMap(m);
}


void Foam::fixedPoroPotentialFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fixedPoroPotentialFvPatchScalarField& tiptf =
        refCast<const fixedPoroPotentialFvPatchScalarField>(ptf);

    h0_.rmap(tiptf.h0_, addr);
}


void Foam::fixedPoroPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const poroHydraulicModel &poroHydraulic_ = this->db().time().lookupObject<poroHydraulicModel>("poroHydraulicModel");
    h0_ = headSeries_->value(this->db().time().timeOutputValue());

        //Info << "h at patch " << patch().name() << " = "
        //    << headSeries_(this->db().time().timeOutputValue())
        //     << endl;
    
if (isHead_)
{
        const fvPatchField<scalar>& z =
            poroHydraulic_.z().boundaryField()[this->patch().index()];
    operator==
    (
        h0_-z
    );
}
else
{
    operator==
    (
        (h0_ - poroHydraulic_.href().value())*poroHydraulic_.magGamma().value()
    );
}

    fixedValueFvPatchScalarField::updateCoeffs();
}



void Foam::fixedPoroPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    h0_.writeEntry("h0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedPoroPotentialFvPatchScalarField
    );
}

// ************************************************************************* //
