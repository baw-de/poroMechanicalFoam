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
    solids4Foam

Description
    Initialization of stress field, that includes the pore water.
    Flow solver will be initalized to get all saturation.

    Then the solid solver will calculate a consolidation under total soil weight. 

Author
    Denis Maier, BAW.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "poroFluidModel.H"
#include "varSatPoroHydraulicModel.H"
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

    // Create the general physics class
    autoPtr<poroFluidModel> fluid = poroFluidModel::New(runTime,"poroFluid",false);
    
    autoPtr<solidModel> solid = solidModel::New(runTime,"solid");

        #include "createMeshToMeshInterpolation.H"

        const varSatPoroHydraulicModel& poroHydraulic = fluid().lookupObject<varSatPoroHydraulicModel>("poroHydraulicProperties");

        fluid().p() = fluid().p_rgh() + poroHydraulic.p_Hyd();
        forAll(fluid().p().boundaryField(), iPatch)
        {
            fluid().p().boundaryFieldRef()[iPatch] =
                fluid().p_rgh().boundaryField()[iPatch] 
                + poroHydraulic.p_Hyd().boundaryField()[iPatch];
        }
        poroHydraulic.S(fluid().p());   

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
                "p",
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
            solidToPoroFluid_().mapSrcToTgt(poroHydraulic.S())
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
            solidToPoroFluid_().mapSrcToTgt(poroHydraulic.n())
        );


        solid.ref().recalculateRho();

        const solidModel& solidRef(solid());
        volScalarField& rho = const_cast<volScalarField&>(solidRef.rho());
        volScalarField rhoEnd = rho;
        
    while(runTime.run())
    {
        runTime++;

        rho = runTime.timeOutputValue()*rhoEnd;
        Info << "current density: " << max(rho) << endl;

        solid().evolve();

        solid().updateTotalFields();

        solid().writeFields(runTime);
        solidRef.rho().write();
        runTime.writeNow();
    }

    solid().end();
    fluid().end();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
