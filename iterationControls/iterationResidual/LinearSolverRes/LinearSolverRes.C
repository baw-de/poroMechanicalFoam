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

#include "LinearSolverRes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam{
     defineTypeNameAndDebug(LinearSolverRes, 0);

    addToRunTimeSelectionTable(
            iterationResidual,
            LinearSolverRes,
            dictionary);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LinearSolverRes::LinearSolverRes
(
    const Time& runTime,
    const word name,
    const ITstream stream,
    const bool writeField
)
:
    iterationResidual(runTime, name, stream, writeField)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LinearSolverRes::~LinearSolverRes()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::LinearSolverRes::calcResidual()
{
    List<scalar> linSolvResiduals;
    HashTable<const fvMesh*> meshes = runTime().lookupClass<const fvMesh>();
    forAllConstIters(meshes, meshIter)
    {
        const fvMesh& regionMesh = *(meshIter.val());
            forAllConstIters(regionMesh.solverPerformanceDict(), entryIter)
            {
            List<Foam::SolverPerformance<double>> sp(entryIter->stream());
            linSolvResiduals.append(sp.first().initialResidual());
            }
    }

    residual_ = operation(linSolvResiduals);
    return residual_;
}


// * * * * * * * * * * * * * * Ostream operation  * * * * * * * * * * * * * * //

// ************************************************************************* //
