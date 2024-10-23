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
    poroMechanicalFoam

Description
    specialist solver for geotechnical problems, based on solids4foam,
    where the solved mathematical model (fluid, solid or fluid-solid)
    is chosen at run-time.

    We choose a physicsModel which might either be
    solid -> pure solid, e.g. dry / drained
    poroFluid -> pure groundwater 
    poroSolid -> coupled flow deformation

    Depending on the choice the physicsModel object is 
    INTERNALLY (not known to this application, the application
    only knows the interfaces provided by the physicsModel class)
    one of the following:
    -> solidModel
    -> poroFluidModel
    -> poroSolidInterface

    We call the "new" function of one of these classes.

    We then proceed with the time loop. where we periodically call the
    evolve function of the physicsModel object (not class! internally
    this is the run-time chosen model and the physicsModel evolve function 
    is overwritten by the chosen model evolve function)

Author
    Denis Maier, BAW.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "physicsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "listAllRegistered",
        "List all the registered objects"
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "solids4FoamWriteHeader.H"
#   include "postProcess.H"  // makes it possible to use postProcessing tools while running

    // Should we adjust the timestep or is it constant?
    bool adjustTimeStep =
        runTime.controlDict().getOrDefault("adjustTimeStep", false);

    // Create the general physics object that internally knows how to handly requests 
    // depending on the run-time chosen model
    autoPtr<physicsModel> physics = physicsModel::New(runTime);

    // Write all registered fields if we added the command line argument
    // So that we know whats in memory.
    if (args.optionFound("listAllRegistered"))
    {
        // Get a list of all the names in of all registered fields
        wordList regs(runTime.sortedNames());
        // The runTime has a tree structure, so fields might by at multiple levels
        // We first write out the first level
        Info << "Objects in " << runTime.name() << ":" << nl
             <<  regs << endl;
        
        // We loop though the entries of the first level
        forAll(regs,iReg)
        {
            // if the entry is a registry (branch of the tree) itself 
            // We write out all the entries of this branch 
            // Otherwise the entry was a field and we skip to the next entry
            if (runTime.foundObject<objectRegistry>(regs[iReg]))
            {
                Info << "Objects in " << regs[iReg] << ":" << nl
                     << runTime.subRegistry(regs[iReg]).sortedNames() << endl;
            }
        }
    }

    ///-------  Starting the time loop!  ---------///
    while (runTime.run())
    {
        // Update deltaT, if desired, before moving to the next step
        if (adjustTimeStep){
            physics().setDeltaT(runTime);
        }
        
        // Advancing the time
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Solve the mathematical model
        physics().evolve();

        // Let the physics model know the end of the time-step has been reached
        // To calculate some fields we only want to calculate at the end of the 
        // time step
        physics().updateTotalFields();

        // Check if this is one of the timesteps that we want to write fields to disk
        if (runTime.outputTime())
        {
            // Write all registered and AUTO_WRITE marked fields to disk
            physics().writeFields(runTime);
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    ///-------  End of time loop!  ---------///
    physics().end();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
