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

#include "brooksCorey.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace saturationLaws
  {
    defineTypeNameAndDebug(brooksCorey, 0);

    addToRunTimeSelectionTable(
        saturationLaw,
        brooksCorey,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
    scalar brooksCorey::H
    (
      const scalar p,
      const scalar pD,
      const scalar n
    ) const
    {
        if (p <= pD)
        {
            return (pow(pD / p, n));
        }
        else
        {
            return 1.0;
        }
    }

    scalar brooksCorey::S
    (
      const scalar p,
      const scalar pD,
      const scalar n,
      const scalar S0,
      const scalar Sr
    ) const
    {
        return H(p,pD,n) * (S0 - Sr) + Sr;
    }

    scalar brooksCorey::krFunc
    (
      const scalar p,
      const scalar pD,
      const scalar n
    ) const
    {
      return pow(H(p,pD,n), 3 + 2 / n);
    }

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    brooksCorey::brooksCorey(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : saturationLaw(name, poroHydraulicProperties, pField),
          brooksCoreyCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),

          n_(
              IOobject(
                  "n_swcc",
                  pField.time().timeName(),
                  db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              dimensionedScalar(brooksCoreyCoeffs_.lookup("n"))),

          S_r(
              IOobject(
                  "S_r",
                  pField.time().timeName(),
                  db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              dimensionedScalar(brooksCoreyCoeffs_.lookup("S_r"))),
          S_0(
              IOobject(
                  "S_pe",
                  pField.time().timeName(),
                  db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              dimensionedScalar(brooksCoreyCoeffs_.lookup("S_pe"))),
          pD_(
              IOobject(
                  "p_e",
                  pField.time().timeName(),
                  db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              dimensionedScalar(brooksCoreyCoeffs_.lookup("p_e")))

    {
      if (pD_.dimensions()!=pField.dimensions())
      {
        FatalErrorIn("vanGenuchten::vanGenuchten") << "Pressure dimensions are: " << pField.dimensions() << endl
                                                   << "p_e has wrong dimensions " << pD_.dimensions() << endl;
      }
      if (debug)
      {
        pD_.write();
        S_0.write();
        S_r.write();
        n_.write();
      }
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar brooksCorey::pStar(const label cellI) const
    {
      return pD_.internalField()[cellI];
    }

    scalar brooksCorey::C(const scalar p, const label cellI)
    {
      scalar pD = pD_.internalField()[cellI];
        if (p <= pD)
        {
          scalar n = n_.internalField()[cellI];
          scalar S0 = S_0.internalField()[cellI];
          scalar Sr = S_r.internalField()[cellI];
          return pow(pD / p, n + 1) * n * (S0 - Sr) / mag(pD);
        }
        else
        {
          return 0.0;
        }
    }

    scalar brooksCorey::S(const scalar p, const label cellI)
    {
        scalar pD = pD_.internalField()[cellI];
        scalar n = n_.internalField()[cellI];
        scalar S0 = S_0.internalField()[cellI];
        scalar Sr = S_r.internalField()[cellI];

        return S(p,pD,n,S0,Sr);
    }

    scalar brooksCorey::S
    (
        const scalar p,
        const label patchI,
        const label faceI
    )
    {
        scalar pD = pD_.boundaryField()[patchI][faceI];
        scalar n = n_.boundaryField()[patchI][faceI];
        scalar S0 = S_0.boundaryField()[patchI][faceI];
        scalar Sr = S_r.boundaryField()[patchI][faceI];

        return S(p,pD,n,S0,Sr);
    }

    scalar brooksCorey::kr(const scalar p, const label cellI)
    {
      scalar pD = pD_.internalField()[cellI];
      scalar n = n_.internalField()[cellI];
      return krFunc(p,pD,n);
    }

    scalar brooksCorey::kr
    (
        const scalar p,
        const label patchI,
        const label faceI
    )
    {
      scalar pD = pD_.boundaryField()[patchI][faceI];
      scalar n = n_.boundaryField()[patchI][faceI];
      return krFunc(p,pD,n);
    }


  } // End of namespace saturationLaws
} // End of namespace Foam

//*********************************************************** //
