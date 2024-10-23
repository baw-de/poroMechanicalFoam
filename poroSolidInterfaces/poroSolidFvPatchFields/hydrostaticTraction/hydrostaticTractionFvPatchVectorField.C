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

#include "hydrostaticTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hydrostaticTractionFvPatchVectorField::
hydrostaticTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    poroTractionFvPatchVectorField(p, iF),
    pressure_(p.size(), 0.0),
    hw_(0.0),
    gamma_(0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


hydrostaticTractionFvPatchVectorField::
hydrostaticTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    poroTractionFvPatchVectorField(p, iF, dict),
    pressure_(p.size(), 0.0),
    hw_(readScalar(dict.lookup("hw"))),
    gamma_(dict.lookupOrDefault<scalar>("gamma",9810))
{
    Info<< "Creating " << type() << " boundary condition" << endl;

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

}


hydrostaticTractionFvPatchVectorField::
hydrostaticTractionFvPatchVectorField
(
    const hydrostaticTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    poroTractionFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    pressure_(mapper(stpvf.pressure_)),
#else
    pressure_(stpvf.pressure_, mapper),
#endif
    hw_(stpvf.hw_),
    gamma_(stpvf.gamma_)
{}


hydrostaticTractionFvPatchVectorField::
hydrostaticTractionFvPatchVectorField
(
    const hydrostaticTractionFvPatchVectorField& stpvf
)
:
    poroTractionFvPatchVectorField(stpvf),
    pressure_(stpvf.pressure_),
    hw_(stpvf.hw_),
    gamma_(stpvf.gamma_)
{}


hydrostaticTractionFvPatchVectorField::
hydrostaticTractionFvPatchVectorField
(
    const hydrostaticTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    poroTractionFvPatchVectorField(stpvf, iF),
    pressure_(stpvf.pressure_),
    hw_(stpvf.hw_),
    gamma_(stpvf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrostaticTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    poroTractionFvPatchVectorField::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(pressure_, pressure_);
#else
    pressure_.autoMap(m);

#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void hydrostaticTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    poroTractionFvPatchVectorField::rmap(ptf, addr);

    const hydrostaticTractionFvPatchVectorField& dmptf =
        refCast<const hydrostaticTractionFvPatchVectorField>(ptf);

    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void hydrostaticTractionFvPatchVectorField::updateCoeffs()
{
    scalarField z(patch().size(),0.0);
    if (updated())
    {
        return;
    }
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    vectorField traction(patch().size(),vector::zero);
    if(solMod.name() == "sixDoFRigidSolid")
    {
        const fvsPatchField<vector>& C = 
        lookupPatchField<surfaceVectorField, vector>("faceLocation");
        z = C.component(vector::Z);
    }
    else
    {
        const vectorField& C = patch().Cf();
        z = C.component(vector::Z);
    }
    
    
    //when mesh is moving this is gonna be implemented instead: scalarField z(patch().Cf().component(vector::Z));
    forAll(patch(),iFace)
    {
        if(z[iFace]<hw_)
        {
            pressure_[iFace] = gamma_ * (hw_ - z[iFace]);
        }
        else
        {
            pressure_[iFace] = 0;
        }
    }

    poroTractionFvPatchVectorField::pressure() = pressure_;

    poroTractionFvPatchVectorField::updateCoeffs();
}


void hydrostaticTractionFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    poroTractionFvPatchVectorField::evaluate();
}


void hydrostaticTractionFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("hw")
        << hw_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma")
        << gamma_ << token::END_STATEMENT << nl;
    poroTractionFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, hydrostaticTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
