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

#include "poroIncrements.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroIncrements, 0);

    addToRunTimeSelectionTable(
        functionObject,
        poroIncrements,
        dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::poroIncrements::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volSymmTensorField &sigma =
            sMesh_().lookupObject<volSymmTensorField>("sigma");

        volSymmTensorField deltaSigma(
            "DSigma",
            sigma - sigma.oldTime());

        deltaSigma.write();

        const volTensorField &gradDD_ =
            sMesh_().lookupObject<volTensorField>("grad(DD)");

        const volSymmTensorField DEpsilon(
            "DEpsilon",
            symm(gradDD_));

        DEpsilon.write();

        const volScalarField &p =
            fMesh_().lookupObject<volScalarField>("p_rgh");

        const volScalarField Dp =
            ("Dp",
             p - p.oldTime());

        Dp.write();
    }

    return true;
}

void Foam::poroIncrements::makeMesh()
{
    word region(dict_.lookupOrDefault<word>("region", "region0"));
    if(runTime_.foundObject<fvMesh>(region))
    {
        fvMesh& mesh(const_cast<fvMesh&>(
            runTime_.lookupObject<fvMesh>(
                 region
                )
            ));
        sMesh_.reset(mesh);
        fMesh_.reset(mesh);
    }
    else
    {
        makeSMesh();
        makeFMesh();
    }
}

void Foam::poroIncrements::makeSMesh()
{
    sMesh_.reset(
                const_cast<fvMesh&>(
                    runTime_.lookupObject<fvMesh>(
                        dict_.lookupOrDefault<word>("poroSolidRegion", "solid")
                    )
                )
            );
}

void Foam::poroIncrements::makeFMesh()
{
    fMesh_.reset(
                const_cast<fvMesh&>(
                    runTime_.lookupObject<fvMesh>(
                     dict_.lookupOrDefault<word>("poroFluidRegion", "poroFluid")
                    )
                )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poroIncrements::poroIncrements(
    const word &name,
    const Time &t,
    const dictionary &dict)
    : functionObject(name),
      name_(name),
      dict_(dict),
      runTime_(t),
      sMesh_(),
      fMesh_()
{
    makeMesh();
    Info << "Creating " << this->name() << " function object" << nl
         << "Using region " << sMesh_().name() << "for poroSolid region" << nl
         << "Using region " << fMesh_().name() << "for poroFluid region" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::poroIncrements::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}

bool Foam::poroIncrements::execute()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}

bool Foam::poroIncrements::read(const dictionary &dict)
{
    return true;
}

bool Foam::poroIncrements::write()
{
    return writeData();
}

// ************************************************************************* //
