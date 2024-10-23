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

#include "limitedGradient.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace conductivityModels
  {
    defineTypeNameAndDebug(limitedGradient, 0);

    addToRunTimeSelectionTable(
        conductivityModel,
        limitedGradient,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    scalar limitedGradient::kFunc(const scalar k0,const scalar iz) const
    {
      return min
      (
        k0 + pos(iz - iCrit_.value()) * kappa_.value() * (1.0 - iCrit_.value()/(iz+SMALL)),
        10*k0
      ); 
    }

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    limitedGradient::limitedGradient(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : conductivityModel(name, poroHydraulicProperties, pField),
        liquefyingZonePermeabilityProperties(poroHydraulicProperties.subDict(typeName + "Coeffs")),
        kappa_(liquefyingZonePermeabilityProperties.lookup("kappa")),
        iCrit_(liquefyingZonePermeabilityProperties.lookup("critGradient")),
        i_(db().objectRegistry::lookupObject<volVectorField>("i")),
        gamma_(db().objectRegistry::lookupObject<uniformDimensionedVectorField>("gamma_water")),
        k0_
        ( 
          IOobject
          (
            "k",
            i_.time().timeName(),
            db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
          ),
          i_.mesh(),
          poroHydraulicProperties.get<dimensionedScalar>("k")      
        )
    {
      if(iCrit_.dimensions() != i_.dimensions() || kappa_.dimensions() != k0_.dimensions())
      {
        FatalErrorInFunction() << "limitedGradient coeffs don't have the right dimensions!" << endl;
      }
      k0_.rename("k0");
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
        
    scalar limitedGradient::k(const label cellI)
    {
      scalar i_z = (i_.internalField()[cellI] & vector(gamma_.value()).normalise());
      return kFunc(k0_.primitiveField()[cellI],i_z);
    }

    scalar limitedGradient::k
    (
        const label patchI,
        const label faceI
    )
    {
      scalar i_z = (i_.boundaryField()[patchI][faceI] & vector(gamma_.value()).normalise());
      return kFunc(k0_.boundaryField()[patchI][faceI],i_z);
    }

  } // End of namespace conductivityModels
} // End of namespace Foam

//*********************************************************** //
