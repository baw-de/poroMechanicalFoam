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

#include "stressChange.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stressChange, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        stressChange,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::stressChange::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volSymmTensorField& sigma =
            mesh_.lookupObject<volSymmTensorField>("sigma");

            volSymmTensorField deltaSigma
                  (
                     "deltaSigma",
                     sigma - sigma0_
                  );

                  deltaSigma.write();

        const volVectorField& D =
            mesh_.lookupObject<volVectorField>("D");

            volVectorField deltaD
                  (
                     "deltaD",
                     D - D0_
                  );

                  deltaD.write();
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stressChange::stressChange
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    ),
    sigma0_(
      IOobject
      (
          "sigma",
          "0",
          mesh_,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      ),
      mesh_
    ),
      D0_(
        IOobject
        (
            "D",
            "0",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
      mesh_
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::stressChange::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}



bool Foam::stressChange::execute()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::stressChange::read(const dictionary& dict)
{
    return true;
}


bool Foam::stressChange::write()
{
    return writeData();
}

// ************************************************************************* //
