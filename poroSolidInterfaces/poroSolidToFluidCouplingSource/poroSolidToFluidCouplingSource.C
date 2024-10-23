/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "poroSolidToFluidCouplingSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "Constant.H"
#include "addToRunTimeSelectionTable.H"
#include "poroSolidInterface.H"


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //
namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(poroSolidToFluidCouplingSource, 0);
        addToRunTimeSelectionTable(option, poroSolidToFluidCouplingSource, dictionary);
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::poroSolidToFluidCouplingSource::poroSolidToFluidCouplingSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    interfaceName_("poroCouplingProperties")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::poroSolidToFluidCouplingSource::addSup
(
        fvMatrix<scalar>& eqn,
        const label fieldi
)
{
    this->addSup(volScalarField::null(), eqn, fieldi);
}


void Foam::fv::poroSolidToFluidCouplingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    // Get poroSold base class, so we can get the coupling terms from it
    const auto& interface = mesh_.thisDb().parent().lookupObject<poroSolidInterface>(interfaceName_,true);
    // we need 1/deltaT because the source should be per unit time. the fixed stress coefficient is not per unit time.
    auto rDeltaT = 1.0/mesh_.time().deltaT();

    // Fixed Stress coefficient, added to the diagonal (Sp part)
    // and RHS with last iteration values(Su part)
    eqn += fvm::SuSp(rDeltaT*interface.implicitCouplingDtoP(),interface.pField());               

    // Explicit coupling (doesnt depend on p, so no implicit parts possible)
    eqn -= fvm::Su(interface.explicitCouplingDtoP(),interface.pField());
}

void Foam::fv::poroSolidToFluidCouplingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
)
{
    this->addSup(volScalarField::null(), eqn, label(1));
}       


bool Foam::fv::poroSolidToFluidCouplingSource::read(const dictionary& dict)
{
    const auto& interface = mesh_.thisDb().parent().lookupObject<poroSolidInterface>(interfaceName_,true);

    fieldNames_.resize(1, interface.pField().name());

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    return true;
}


// ************************************************************************* //
