/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "porosityFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace conductivityModels
  {
    defineTypeNameAndDebug(porosityFunction, 0);

    addToRunTimeSelectionTable(
        conductivityModel,
        porosityFunction,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    porosityFunction::porosityFunction(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : conductivityModel(name, poroHydraulicProperties, pField),
          porosityFunctionProperties(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          kNFunc_(Function1<scalar>::New("kOfN",porosityFunctionProperties))
    {}

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar porosityFunction::k(const label cellI)
    {
      const volScalarField& n_ = db().objectRegistry::lookupObject<volScalarField>("n");
      return kNFunc_().value(n_.internalField()[cellI]);
    }

    scalar porosityFunction::k
    (
        const label patchI,
        const label faceI
    )
    {
      const volScalarField& n_ = db().objectRegistry::lookupObject<volScalarField>("n");
      return kNFunc_().value(n_.boundaryField()[patchI][faceI]);
    }

  } // End of namespace conductivityModels
} // End of namespace Foam

//*********************************************************** //
