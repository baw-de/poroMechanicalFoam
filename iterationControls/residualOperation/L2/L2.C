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
    under the teL2 of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "L2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam
{
    namespace residualOperations
    {
    defineTypeNameAndDebug(L2, 0);

    addToRunTimeSelectionTable(
            residualOperation,
            L2,
            dictionary);
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::residualOperations::L2::L2
(
    const word operation
)
:
    residualOperation(operation)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::residualOperations::L2::~L2()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::residualOperations::L2::operation(const scalarField& x) const
{
    return pow(gSum(pow(x,2)),0.5);
}

Foam::scalar Foam::residualOperations::L2::operation(const List<scalar>& x) const
{
    scalar returnValue;
    forAll(x, ix)
    {
        returnValue += pow(x[ix],2);
    }
    returnValue = pow(returnValue,0.5);
    return returnValue;
}

// ************************************************************************* //
