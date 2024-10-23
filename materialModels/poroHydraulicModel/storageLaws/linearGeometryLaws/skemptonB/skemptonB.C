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

#include "skemptonB.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace storageLaws
  {
    defineTypeNameAndDebug(skemptonB, 0);

    addToRunTimeSelectionTable(
        storageLaw,
        skemptonB,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    skemptonB::skemptonB(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : storageLaw(name, poroHydraulicProperties, pField),
          skemptonBCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          BParam_(
              IOobject(
                  "BParam",
                  mesh().time().timeName(),
                  pField.db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              skemptonBCoeffs_.get<dimensionedScalar>("BParameter")),
          Km_(
              IOobject(
                  "Km",
                  mesh().time().timeName(),
                  pField.db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              skemptonBCoeffs_.get<dimensionedScalar>("Km"))
    {
      if (Km_.dimensions()!=pField.dimensions())
      {
        FatalErrorIn("skemptonB::skemptonB") << "Pressure dimensions are: " << pField.dimensions() << endl
                                             << "Km has wrong dimensions " << Km_.dimensions()  << endl;
      }
      makeCs();
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar skemptonB::Ss
    (
        const scalar n,
        const scalar p,
        const label cellI
    )
    {
      return 
      (
        (
          (1.0 - BParam_.internalField()[cellI]) 
          / (BParam_.internalField()[cellI] * Km_.internalField()[cellI])
        ) 
        + (1.0-n)*Cs_().internalField()[cellI]
      );
    }
  } // End of namespace storageLaws
} // End of namespace Foam

//*********************************************************** //
