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

#include "solidInvariants.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidInvariants, 0);

    addToRunTimeSelectionTable(
        functionObject,
        solidInvariants,
        dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidInvariants::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volSymmTensorField &sigma =
            sMesh_().lookupObject<volSymmTensorField>("sigma");

        const volScalarField P(
            IOobject::groupName("sigma","P"),
            1/3*tr(sigma));

        P.write();

        const volScalarField Q(
            IOobject::groupName("sigma","Q"),
            Foam::sqrt(3.0/2.0*dev(sigma)&&dev(sigma)));

        Q.write();

        const volTensorField &gradDD =
            sMesh_().lookupObject<volTensorField>("grad(DD)");

        const volSymmTensorField symmGradDD(symm(gradDD));

        const volScalarField epsilonDotP(
            IOobject::groupName("DEpsilon","P"),
            1/3*tr(symmGradDD));

        epsilonDotP.write();

        const volScalarField epsilonDotQ(
            IOobject::groupName("DEpsilon","Q"),
            Foam::sqrt(3.0/2.0*dev(symmGradDD)&&dev(symmGradDD)));
        epsilonDotQ.write();
    }

    return true;
}

void Foam::solidInvariants::makeMesh()
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
    }
    else
    {
        makeSMesh();
    }
}

void Foam::solidInvariants::makeSMesh()
{
    sMesh_.reset(
                const_cast<fvMesh&>(
                    runTime_.lookupObject<fvMesh>(
                        dict_.lookupOrDefault<word>("poroSolidRegion", "solid")
                    )
                )
            );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidInvariants::solidInvariants(
    const word &name,
    const Time &t,
    const dictionary &dict)
    : functionObject(name),
      name_(name),
      dict_(dict),
      runTime_(t),
      sMesh_()
{
    makeMesh();
    Info << "Creating " << this->name() << " function object" << nl
         << "Using region " << sMesh_().name() << "for poroSolid region" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidInvariants::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}

bool Foam::solidInvariants::execute()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}

bool Foam::solidInvariants::read(const dictionary &dict)
{
    return true;
}

bool Foam::solidInvariants::write()
{
    return writeData();
}

// ************************************************************************* //
