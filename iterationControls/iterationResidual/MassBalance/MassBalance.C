/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MassBalance.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam{
     defineTypeNameAndDebug(MassBalance, 0);

    addToRunTimeSelectionTable(
            iterationResidual,
            MassBalance,
            dictionary);
}

void Foam::MassBalance::makeMassBalanceRef()
{
    HashTable<const poroFluidModel*> models = runTime().lookupClass<poroFluidModel>();
    bool found = false;
        forAllConstIters(models, modelIter)
        {
            const fvMesh& regionMesh = *(modelIter.val().mesh());
            if(regionMesh.objectRegistry::foundObject<volScalarField>("MassBalanceResidual"))
            {
                found = true;
                MassBalance_.reset(regionMesh.objectRegistry::lookupObject<volScalarField>("MassBalanceResidual"));
            }
        }
    if(!found)
    {
        FatalErrorInFunction << "MassBalance field not found!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MassBalance::MassBalance
(
    const Time& runTime,
    const word name,
    const ITstream stream,
    const bool writeField
)
:
    iterationResidual(runTime, name, stream, writeField),
    writeField_(writeField),
    MassBalance_(),
    HashTable<const poroFluidModel*> models = runTime().lookupClass<poroFluidModel>()
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MassBalance::~MassBalance()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::MassBalance::calcResidual()
{
    if(!MassBalance_.valid())
    {
        makeMassBalanceRef();
    }

    scalarField MassBalanceInternal = mag(MassBalance_().primitiveField());
    if(relative_)
    {
        dimensionedScalar dimensionedSmall("",MassBalance_().dimensions(),SMALL);
        MassBalanceInternal =  mag(MassBalanceInternal)/max(MassBalance_().oldTime(),dimensionedSmall)().primitiveField();
    }
    else
    {

    }
    residual_ = operation(MassBalanceInternal);
    return residual_;
}

void Foam::MassBalance::reset()
    {
        if(!MassBalance_.valid())
        {
            makeMassBalanceRef();
        }
        const_cast<volScalarField&>(MassBalance_()).storeOldTime();
        residual_ = GREAT;
    }

// * * * * * * * * * * * * * * Ostream operation  * * * * * * * * * * * * * * //

// ************************************************************************* //
