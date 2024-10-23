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
    set geostatic inital conditions with a k0-condition.

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // define allowable options and arguments
    argList::validArgs.append("gamma");
    argList::validArgs.append("K0");
    argList args(argc, argv);

#include "createTime.H"
#include "createMesh.H"

    const scalar gamma(readScalar(IStringStream(args.args()[1])()));
    const scalar K0(readScalar(IStringStream(args.args()[2])()));
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const uniformDimensionedVectorField g // gravity
        (
            IOobject(
                "g",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));
    /*if (vector(0, 0, -1) != (g / max(Foam::mag(g), dimensionedScalar("", dimVelocity / dimTime, VSMALL)).value())))
        {
            FoamWarning() << "Mesh axes need to be principle axes and z must be the direction of gravity!!!!" << endl;
        }*/

    volScalarField z = -(mesh.C() & (g / max(Foam::mag(g), dimensionedScalar("", dimVelocity / dimTime, VSMALL))));

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    volSymmTensorField sigma // stress
        (
            IOobject(
                "sigma",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedSymmTensor("0", dimPressure, symmTensor::zero));

    forAll(sigma.internalField(), icell)
    {
        sigma.internalField()[icell] = symmTensor(K0 * gamma * z.internalField()[icell], 0, 0, K0 * gamma * z.internalField()[icell], 0, gamma * z.internalField()[icell]);
    }
    forAll(sigma.boundaryField(), ipatch)
    {
        symmTensorField &sigmaP = sigma.boundaryField()[ipatch];
        scalarField &zP = z.boundaryField()[ipatch];
        forAll(sigmaP, iface)
        {
            sigmaP[iface] = symmTensor(K0 * gamma * zP[iface], 0, 0, K0 * gamma * zP[iface], 0, gamma * zP[iface]);
        }
    }

    sigma -= max(sigma);

    sigma.write();

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
