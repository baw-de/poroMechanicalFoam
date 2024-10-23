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
    helper utitlity to generate a new field empty field in the 0 directory.

Author
    Denis Maier, BAW.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
argList::validArgs.append("kind");
argList::validArgs.append("field");

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    const word kind(word(IStringStream(args.args()[1])()));
    const word field(word(IStringStream(args.args()[2])()));

if (kind == "scalar")
{
    volScalarField vf_(
        IOobject(
            field,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
	dimensionedScalar("",dimless,0.0));
    Info << "Writing " << kind << "Field " << field << " to "
         << vf_.path() << endl;
    vf_.write();
}
else if (kind == "vector")
{
    volVectorField vf_(
        IOobject(
            field,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
	dimensionedVector("",dimless,vector::zero));
    vf_.write();
}
else if (kind == "symmTensor")
{
    volSymmTensorField vf_(
        IOobject(
            field,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
	dimensionedSymmTensor("",dimless,symmTensor::zero));
    vf_.write();
}
else if (kind == "tensor")
{
    volTensorField vf_(
        IOobject(
            field,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
	dimensionedTensor("",dimless,tensor::zero));
    vf_.write();
}
else
{
	FatalError("") << "first argument (kind) must be scalar, vector, symmTensor or tensor"
	<< endl;
}
    Info << nl << "End" << nl << endl;

    return (0);
}

// ************************************************************************* //
