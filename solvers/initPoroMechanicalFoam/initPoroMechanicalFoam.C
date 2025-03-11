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

Application
    initPoroMechanicalFoam

Description
    Initialization of stress field, that includes the pore water.
    Flow solver will be initialized to get all saturation.

    Then the solid solver will calculate a consolidation under total soil weight.

Author
    Denis Maier, BAW.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "poroFluidModel.H"
#include "solidModel.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "gravityLoading",
        "ramp the gravity up slowly"
    );

#   include "setRootCase.H"
#   include "createTime.H"

    // Create the general physics classes
    autoPtr<poroFluidModel> fluid = poroFluidModel::New(runTime, "poroFluid", false);
    autoPtr<solidModel> solid = solidModel::New(runTime, "solid");

#   include "createMeshToMeshInterpolation.H"

    // Get reference to poroPhydraulic model
    const Foam::poroHydraulicModel& poroHydraulic =
        fluid().mesh().lookupObject<poroHydraulicModel>("poroHydraulicModel");

    // Initialize pressure
    fluid().p() = fluid().p_rgh() + poroHydraulic.p_Hyd();
    forAll(fluid().p().boundaryField(), iPatch)
    {
        fluid().p().boundaryFieldRef()[iPatch] =
            fluid().p_rgh().boundaryField()[iPatch]
            + poroHydraulic.p_Hyd().boundaryField()[iPatch];
    }

    // Get saturation
    tmp<volScalarField> saturationField;
    if (fluid().mesh().objectRegistry::foundObject<volScalarField>("S"))
    {
        // If S exists, create a tmp from a reference
        const volScalarField& existingS =
            fluid().mesh().lookupObject<volScalarField>("S");
        saturationField = tmp<volScalarField>(new volScalarField(existingS));
    }
    else
    {
        // Create a new field with default value of 1.0
        saturationField = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "S",
                    runTime.timeName(),
                    fluid().mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid().mesh(),
                dimensionedScalar("S", dimless, 1.0)
            )
        );
    }

    // Map fields from fluid mesh to solid mesh
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            solid().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidToPoroFluid_().mapSrcToTgt(fluid().p())
    );

    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            solid().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidToPoroFluid_().mapSrcToTgt(fluid().p_rgh())
    );

    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            solid().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidToPoroFluid_().mapSrcToTgt(saturationField())
    );

    volScalarField n
    (
        IOobject
        (
            "n",
            runTime.timeName(),
            solid().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidToPoroFluid_().mapSrcToTgt(poroHydraulic.n0())
    );

    // Recalculate density for the solid model
    solid.ref().recalculateRho();

    const solidModel& solidRef(solid());
    volScalarField& rho = const_cast<volScalarField&>(solidRef.rho());
    volScalarField rhoEnd = rho;

    while(runTime.run())
    {
        runTime++;

        // Ramp up the density over time
        rho = runTime.timeOutputValue() * rhoEnd;
        Info << "Current density: " << max(rho) << endl;

        solid().evolve();
        solid().updateTotalFields();

        // Write out fields
        solid().writeFields(runTime);
        solidRef.rho().write();
        runTime.writeNow();
    }

    solid().end();
    fluid().end();

    Info<< nl << "End" << nl << endl;

    return(0);
}
