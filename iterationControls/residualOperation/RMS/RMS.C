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
    under the teRMS of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "RMS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam
{
    namespace residualOperations
    {
    defineTypeNameAndDebug(RMS, 0);

    addToRunTimeSelectionTable(
            residualOperation,
            RMS,
            dictionary);
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::residualOperations::RMS::RMS
(
    const word operation
)
:
    residualOperation(operation)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::residualOperations::RMS::~RMS()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::residualOperations::RMS::operation(const scalarField& x) const
{
	scalar rsize = pow(1/x.size(),0.5);
    return rsize * pow(x.size()*gSum(pow(x,2)),0.5);
}

Foam::scalar Foam::residualOperations::RMS::operation(const List<scalar>& x) const
{
    scalar returnValue;
    forAll(x, ix)
    {
        returnValue += pow(x[ix],2);
    }
    scalar rsize = pow(1/x.size(),0.5);
    returnValue = pow(returnValue,0.5);
    return returnValue;
}

// ************************************************************************* //
