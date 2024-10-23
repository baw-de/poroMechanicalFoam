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

#include "kozenyCarman.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace conductivityModels
  {
    defineTypeNameAndDebug(kozenyCarman, 0);

    addToRunTimeSelectionTable(
        conductivityModel,
        kozenyCarman,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
    scalar kozenyCarman::kFunc(const scalar nScalar) const
    {
      return 
      (
        pow(nScalar, 3) * pow(D50_.value(), 2) / (mu_.value() * 180.0 * pow(1 - nScalar, 2))
      );
    }
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    kozenyCarman::kozenyCarman(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : conductivityModel(name, poroHydraulicProperties, pField),
          kozenyCarmanProperties(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          D50_(kozenyCarmanProperties.lookup("D50")),
          mu_(kozenyCarmanProperties.lookup("viscosity"))
    {}

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar kozenyCarman::k(const label cellI)
    {
      const volScalarField& nk = db().objectRegistry::lookupObject<volScalarField>("n");
      const scalar n_ = nk.internalField()[cellI];
      return kFunc(n_);
    }

    scalar kozenyCarman::k
    (
        const label patchI,
        const label faceI
    )
    {
      const volScalarField& nk = db().objectRegistry::lookupObject<volScalarField>("n");
      const scalar n_ = nk.boundaryField()[patchI][faceI];
      return kFunc(n_);
    }

  } // End of namespace conductivityModels
} // End of namespace Foam

//*********************************************************** //
