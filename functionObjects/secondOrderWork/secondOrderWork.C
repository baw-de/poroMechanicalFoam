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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "secondOrderWork.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(secondOrderWork, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        secondOrderWork,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::secondOrderWork::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volSymmTensorField& sigma =
            runTime_.lookupObject<volSymmTensorField>("sigma");

            volSymmTensorField deltaSigma
                  (
                     "deltaSigma",
                     sigma - sigma.oldTime()
                  );

                  deltaSigma.write();

        const volTensorField& gradDD_ =
            runTime_.lookupObject<volTensorField>("grad(DD)");

        const volScalarField d2W
                  (
                     "d2W",
                     deltaSigma && symm(gradDD_)
                  );

                  d2W.write();
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::secondOrderWork::secondOrderWork
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t)
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::secondOrderWork::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}



bool Foam::secondOrderWork::execute()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::secondOrderWork::read(const dictionary& dict)
{
    return true;
}


bool Foam::secondOrderWork::write()
{
    return writeData();
}

// ************************************************************************* //
