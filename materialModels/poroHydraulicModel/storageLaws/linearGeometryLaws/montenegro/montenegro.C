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

#include "montenegro.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace storageLaws
  {
    defineTypeNameAndDebug(montenegro, 0);

    addToRunTimeSelectionTable(
        storageLaw,
        montenegro,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    montenegro::montenegro(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : storageLaw(name, poroHydraulicProperties, pField),
          montenegroCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          S_pe_(
              IOobject(
                  "S_pe",
                  mesh().time().timeName(),
                  pField.db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              montenegroCoeffs_.get<dimensionedScalar>("S_e")),
          p_At_(montenegroCoeffs_.get<dimensionedScalar>("p_At")),
          p_e_(montenegroCoeffs_.get<dimensionedScalar>("p_e"))
    {
      if (p_At_.dimensions()!=pField.dimensions() || p_e_.dimensions()!=pField.dimensions())
      {
              FatalErrorIn("montenegro::montenegro") << "Pressure dimensions are: " << pField.dimensions() << endl
                                                    << "one or more coefficients have wrong dimensions " << endl;
      }
      makeCs();
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar montenegro::Ss
    (
        const scalar n,
        const scalar p,
        const label cellI
    )
    {
      if (p > p_e_.value())
      {
        return 
        (
          n * 
          p_At_.value() * (1.0 - S_pe_.primitiveField()[cellI]) 
          / pow(p_At_.value() + p, 2)
          + (1.0-n) * Cs_().primitiveField()[cellI]
        );
      }
      else
      {
        return 0.0;
      }
    }

  } // End of namespace storageLaws
} // End of namespace Foam

//*********************************************************** //
