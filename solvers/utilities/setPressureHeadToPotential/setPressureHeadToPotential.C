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
    set potential (h) by giving pressure head (p/gamma_w).
    This is useful when using the potential / potential pressure solver.

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "UniformDimensionedField.H"
//#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // does it work correctly in parallel? If not you need this:
    // argList::noParallel();

    // define allowable options and arguments
    argList::validOptions.insert("pressureUnits", "Switch");
    argList args(argc, argv);

#include "createTime.H"
#include "createMesh.H"

    Switch inPressureUnits = args.optionLookupOrDefault("pressureUnits", true);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info
        << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    const volScalarField initHead(
        IOobject(
            "pHead",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE),
        mesh);

    //autoPtr<poroHydraulicModel> poroHydraulic(new poroHydraulicModel(initHead));

    IOdictionary poroHydraulicDict(
            IOobject(
                "poroHydraulicProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    dimensionedVector gamma = dimensionedVector(poroHydraulicDict.lookup("gamma"));

    volScalarField z = mesh.C() &  vector(gamma.value()).normalise();

    Info << "Water specific weight has been read to: " << gamma << endl; 


    volScalarField Potential // Pressure head
        (
            IOobject(
                "Potential",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            //poroHydraulic().z() + initHead);
            initHead - z); // use initHead as initializer to copy boundary conditions

    if (!inPressureUnits)
    {
        Warning << "You selected head units for Potential output field!" << endl;
        Potential.write();
    }
    else if (inPressureUnits)
    {
        volScalarField p_rgh // Pressure head
            (
                IOobject(
                    "p_rgh",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                //(Potential - poroHydraulic().href()) * poroHydraulic().magGamma());
                (Potential -  dimensionedScalar(poroHydraulicDict.lookup("href"))) * mag(gamma),
                initHead.boundaryField().types());
        p_rgh.write();
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
