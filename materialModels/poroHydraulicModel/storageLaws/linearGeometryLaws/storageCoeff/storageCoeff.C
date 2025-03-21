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

#include "storageCoeff.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace storageLaws
  {
    defineTypeNameAndDebug(storageCoeff, 0);

    addToRunTimeSelectionTable(
        storageLaw,
        storageCoeff,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    storageCoeff::storageCoeff(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : storageLaw(name, poroHydraulicProperties, pField),
          storageCoeffCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          SsScalar_(storageCoeffCoeffs_.get<dimensionedScalar>("Ss")),
          SsField_()
    {
      IOobject SsHeader(
        "Ss",
        mesh().time().timeName(),
        db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE);

      if(SsHeader.typeHeaderOk<volScalarField>())
      {
        SsField_.reset(
          new volScalarField
          (
            SsHeader,
            mesh())
          );
        if (SsField_().dimensions()!=dimless/pField.dimensions())
        {
                FatalErrorIn("storageCoeff::storageCoeff") << "Pressure dimensions are: " << pField.dimensions() << endl
                                                      << "Ss has wrong dimensions " << SsField_().dimensions()  << endl;
        }
      }
      else
      {
        if (SsScalar_.dimensions()!=dimless/pField.dimensions())
        {
                FatalErrorIn("storageCoeff::storageCoeff") << "Pressure dimensions are: " << pField.dimensions() << endl
                                                      << "Ss has wrong dimensions " << SsScalar_.dimensions()  << endl;
        }
      }

    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    scalar storageCoeff::Ss
    (
        const scalar n,
        const scalar p,
        const label cellI
    )
    {
      if(SsField_.valid())
      {
        return SsField_().primitiveField()[cellI];
      }
      return SsScalar_.value();
    };

  } // End of namespace storageLaws
} // End of namespace Foam

//*********************************************************** //
