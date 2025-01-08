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

\*---------------------------------------------------------------------------*/

#include "[template].H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "iterationControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace poroFluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Register class for runtime selection
defineTypeNameAndDebug([template], 0);
addToRunTimeSelectionTable(poroFluidModel, [template], dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

poroHydraulicModel& [template]::poroHydraulic()
{
    if (!poroHydPtr_.valid())
    {
        FatalErrorIn("[template]::poroHydraulic()")
            << "poroHydraulicModel not initialized before first use" << endl;
    }
    return poroHydPtr_.ref();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

[template]::[template]
(
    Time& runTime,
    const word& region,
    const bool sharedmesh
)
:
    // Initialize base class
    poroFluidModel(typeName, runTime, "p_rgh", region, sharedmesh),
    
    // Create hydraulic model
    poroHydPtr_(new poroHydraulicModel(p_rgh(), g())),
    
    // Initialize storage coefficient field
    Ss_
    (
        IOobject
        (
            "Ss",
            runTime.timeName(),
            p_rgh().db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        poroHydraulic().Ss(n(), p_rgh())
    ),
    
    // Initialize hydraulic conductivity field
    kbyGammaf_
    (
        IOobject
        (
            "kEffbyGammaf",
            runTime.timeName(),
            p_rgh().db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        poroHydraulic().kf() / poroHydraulic().magGamma()
    )
{
    // Ensure p_rgh field exists
    pRGHisRequired();
    
    // Initialize flux field
    phi() = (-kbyGammaf_ * mesh().magSf() * fvc::snGrad(p_rgh()));
    
    Info << "Done generating poroFluidModel" << endl;

    // Write debug fields if requested
    if(debug)
    {
        n().write();
        poroHydraulic().k()().write();
        Ss_.write();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool [template]::evolve()
{
    // Initialize performance tracking
    SolverPerformance<scalar> solverPerfp;
    SolverPerformance<scalar>::debug = 0;
    
    // Reset iteration control
    iterCtrl().reset();

    Info << "Evolving fluid model: " << this->type() << endl;

    do  // Non-orthogonal correctors loop
    {
        // Store previous iteration for relaxation
        p_rgh().storePrevIter();

        if (debug)
        {
            Info << "Assembling pEqn" << endl;
        }

        // Update storage coefficient if needed
        if(poroHydraulic().updatesSs())
        {
            Ss_ = poroHydraulic().Ss(n(), p_rgh());
        }

        // Assemble and solve pressure equation
        fvScalarMatrix pEqn
        (
            fvm::ddt(Ss_, p_rgh())                  // Storage term
          - fvm::laplacian(kbyGammaf_, p_rgh())     // Pressure flux term
         == fvOptions()(Ss_, p_rgh())               // Source terms
        );
        
        // Apply constraints
        fvOptions().constrain(pEqn);

        // Solve system
        solverPerfp = pEqn.solve();

        // Apply relaxation
        p_rgh().relax();
        
        // Apply corrections from fvOptions
        fvOptions().correct(p_rgh());
        
        // Update derived fields
        phi() = pEqn.flux();
        hydraulicGradient() = fvc::grad(p_rgh())/poroHydraulic().magGamma();

        // Update conductivity if needed
        if(poroHydraulic().updatesK())
        {
            kbyGammaf_ = poroHydraulic().kf() / poroHydraulic().magGamma();
        }

        // Update boundary conditions
        p_rgh().correctBoundaryConditions();

    } while (iterCtrl().loop());

    // Write iteration information
    iterCtrl().write();

    Info << "[template] evolved" << endl;

    // Update total pressure and pressure time derivative
    p() = p_rgh() + poroHydraulic().p_Hyd();
    pDot() = fvc::ddt(p());
    
    return 0;
}

scalar [template]::newDeltaT()
{
    // Calculate maximum hydraulic diffusivity
    scalar maxCv = max(poroHydraulic().k() / (Ss_*poroHydraulic().magGamma())).value();
    Info << "max c_v: " << maxCv << endl;
    
    // Calculate minimum cell size
    scalar minDeltaX = min(mesh().V() / fvc::surfaceSum(mag(mesh().Sf())));
    
    // Calculate time step based on Courant number
    scalar wishedTimeStep = CoNumber_ * pow(minDeltaX, 2) / maxCv;
    scalar maxDeltaTFact = wishedTimeStep / runTime().deltaTValue();
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);
    
    return deltaTFact * runTime().deltaTValue();
}

void [template]::writeFields(const Time& runTime)
{
    poroFluidModel::writeFields(runTime);
}

void [template]::end()
{
    poroHydraulic().write();
    poroFluidModel::end();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace poroFluidModels
} // End namespace Foam

// ************************************************************************* //