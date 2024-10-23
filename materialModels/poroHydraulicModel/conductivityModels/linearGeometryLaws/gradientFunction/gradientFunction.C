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

#include "gradientFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace conductivityModels
  {
    defineTypeNameAndDebug(gradientFunction, 0);

    addToRunTimeSelectionTable(
        conductivityModel,
        gradientFunction,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    gradientFunction::gradientFunction(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : conductivityModel(name, poroHydraulicProperties, pField),
          gradientFunctionProperties(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          kIFunc_(Function1<scalar>::New("kOfI",gradientFunctionProperties)),
          i_(db().objectRegistry::lookupObject<volVectorField>("i")),
          gamma_(db().objectRegistry::lookupObject<uniformDimensionedVectorField>("gamma_water"))
    {}

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
    
    scalar gradientFunction::k(const label cellI)
    {
      scalar i_z = (i_.internalField()[cellI] & vector(gamma_.value()).normalise());
      return kIFunc_().value(i_z);
    }

    scalar gradientFunction::k
    (
        const label patchI,
        const label faceI
    )
    {
      scalar i_z = (i_.boundaryField()[patchI][faceI] & vector(gamma_.value()).normalise());
      return kIFunc_().value(i_z);
    }

  } // End of namespace conductivityModels
} // End of namespace Foam

//*********************************************************** //
