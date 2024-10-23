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
    UNFINISHED WORK!

Author
    Denis Maier, BAW.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "physicsModel.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "solids4FoamWriteHeader.H"
    IOdictionary mechanicalDict(
        IOobject(
            "mechanicalProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE));
    PtrList<entry> lawEntries(mechanicalDict.lookup("mechanical"));
    dictionary to_reduce_Dict = lawEntries[0].dict();
    dimensionedScalar startPhi_(to_reduce_Dict.lookup("startFrictionAngle"));
    Info << "Starting with Friction angle: " << startPhi_ << endl;
    startPhi_.value() = degToRad(startPhi_.value());
    dimensionedScalar startC_(to_reduce_Dict.lookup("startCohesion"));
    Info << "Starting with cohesion: " << startC_ << endl;
    scalar FoS_(runTime.timeOutputValue() + runTime.deltaTValue());
    // Create the general physics class
    dimensionedScalar newC("cohesion", startC_.dimensions(), startC_.value() / FoS_);
    dimensionedScalar newPhi("frictionAngle", startPhi_.dimensions(), radToDeg(Foam::atan(Foam::tan(startPhi_.value()) / FoS_)));
    dimensionedScalar newPsi("dilationAngle", newPhi);
    to_reduce_Dict.set("frictionAngle", newPhi);
    to_reduce_Dict.set("cohesion", newC);
    to_reduce_Dict.set("dilationAngle", newPsi);
    lawEntries[0].dict() = to_reduce_Dict;
    mechanicalDict.set("mechanical", lawEntries);
    mechanicalDict.regIOobject::write();

    do
    {

        Info << "Current factor of safety: " << FoS_ << endl;
        Info << "Current cohesion: " << newC << endl;
        Info << "Current friction angle: " << newPhi << endl;
        word prevTimeName(runTime.timeName());
        autoPtr<physicsModel> physics = physicsModel::New(runTime);
        runTime++;
        Info << "Current Factor of Safety: " << runTime.timeOutputValue() << endl;
        // Solve the mathematical model
        physics().evolve();
        // Let the physics model know the end of the time-step has been reached
        physics().updateTotalFields();

        physics().writeFields(runTime);
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        physics().end();

        FoS_ = runTime.timeOutputValue() + runTime.deltaTValue();

        newC.value() = startC_.value() / FoS_;
        newPhi.value() = radToDeg(Foam::atan(Foam::tan(startPhi_.value()) / FoS_));
        newPsi.value() = newPhi.value();
        to_reduce_Dict.set("frictionAngle", newPhi);
        to_reduce_Dict.set("cohesion", newC);
        to_reduce_Dict.set("dilationAngle", newPsi);
        lawEntries[0].dict() = to_reduce_Dict;
        mechanicalDict.set("mechanical", lawEntries);
        mechanicalDict.regIOobject::write();
        cp(runTime.path() / prevTimeName / "TotalHead", runTime.path() / runTime.timeName() / "TotalHead");
        cp(runTime.path() / prevTimeName / "Sw", runTime.path() / runTime.timeName() / "Sw");
        physics().setDeltaT(runTime);

    } while (!runTime.end());

    Info
        << nl << "End" << nl << endl;

    return (0);
}

// ************************************************************************* //
