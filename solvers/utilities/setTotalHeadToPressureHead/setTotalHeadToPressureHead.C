/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    set hydrostatic pressure head (p/gamma_w) by giving total head (h).
    This is useful when using the pressure head based flow solver.

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // does it work correctly in parallel? If not you need this:
    // argList::noParallel();

    // define allowable options and arguments
    argList::validOptions.insert("pressureUnits", "Switch");
    argList::validOptions.insert("rho", "scalar");
    argList args(argc, argv);

#include "createTime.H"
#include "createMesh.H"

    Switch inPressureUnits = args.optionLookupOrDefault("pressureUnits", false);
    scalar rho_scalar = args.optionLookupOrDefault("rho", 1000);
    dimensionedScalar rho("rho", dimDensity, rho_scalar);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info
        << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    const uniformDimensionedVectorField g // Pressure head
        (
            IOobject(
                "g",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));

    const volScalarField initHead(
        IOobject(
            "initialTotalHead",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE),
        mesh);

    volScalarField z = (mesh.C() & (g / max(Foam::mag(g), dimensionedScalar("", dimVelocity / dimTime, VSMALL))));

    if (!inPressureUnits)
    {

        const volScalarField pHeadInit // Pressure head
            (
                IOobject(
                    "pHeadTemplate",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE),
                mesh);
        volScalarField pHead // Pressure head
            (
                IOobject(
                    "pHead",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                pHeadInit);

        pHead = z + initHead;

        forAll(pHead.boundaryField(), iPatch)
        {
            forAll(pHead.boundaryField()[iPatch], iFace)
            {
                pHead.boundaryField()[iPatch][iFace] = z.boundaryField()[iPatch][iFace] + initHead.boundaryField()[iPatch][iFace];
            }
        }
        pHead.write();
    }
    else if (inPressureUnits)
    {
        Warning << "You selected pressure units for output field!" << endl;
        const volScalarField pHeadInit // Pressure head
            (
                IOobject(
                    "pTemplate",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE),
                mesh);
        volScalarField pHead // Pressure head
            (
                IOobject(
                    "p",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                pHeadInit);

        pHead = (z + initHead) * (rho * mag(g));

        forAll(pHead.boundaryField(), iPatch)
        {
            forAll(pHead.boundaryField()[iPatch], iFace)
            {
                pHead.boundaryField()[iPatch][iFace] = (z.boundaryField()[iPatch][iFace] + initHead.boundaryField()[iPatch][iFace]) * (rho.value() * mag(g).value());
            }
        }
        pHead.write();
    }
    else
    {
        SeriousErrorIn("setTotalHeadPressueHead") << "Neither PressureHead nor Pressure Units were selected. What happend?" << endl;
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
