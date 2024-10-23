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

#include "KPrime.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace storageLaws
  {
    defineTypeNameAndDebug(KPrime, 0);

    addToRunTimeSelectionTable(
        storageLaw,
        KPrime,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    KPrime::KPrime(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : storageLaw(name, poroHydraulicProperties, pField),
          KPrimeCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          S_p0_(
              IOobject(
                  "S_p0",
                  mesh().time().timeName(),
                  pField.db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              KPrimeCoeffs_.get<dimensionedScalar>("Sw_0")),
          Kw_(KPrimeCoeffs_.get<dimensionedScalar>("Kw")),
          p_At_(KPrimeCoeffs_.get<dimensionedScalar>("p_At")),
          pDep_(Switch(KPrimeCoeffs_.lookupOrAddDefault("pDependent",false)))
    {
      if (p_At_.dimensions()!=pField.dimensions() || Kw_.dimensions()!=pField.dimensions())
      {
              FatalErrorIn("KPrime::KPrime") << "Pressure dimensions are: " << pField.dimensions() << endl
                                                    << "one or more coefficients have wrong dimensions " << endl;
      }
      makeCs();
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar KPrime::Ss
    (
        const scalar n,
        const scalar p,
        const label cellI
    )
    {
        scalar tp = p;
        if (pDep_ != true)
        {
          tp = 0.0;
        }

        return 
        (
          n*
          (
            S_p0_.internalField()[cellI] / Kw_.value() 
            + (1.0 - S_p0_.internalField()[cellI]) / (p_At_.value() + tp)
          ) 
          + (1.0 - n) * Cs_().internalField()[cellI]
        );
    }

  } // End of namespace storageLaws
} // End of namespace Foam
//*********************************************************** //
