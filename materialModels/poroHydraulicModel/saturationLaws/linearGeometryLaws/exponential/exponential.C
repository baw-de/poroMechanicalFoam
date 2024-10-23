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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {
        defineTypeNameAndDebug(exponential, 0);

        addToRunTimeSelectionTable(
            saturationLaw,
            exponential,
            dictionary);

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
        scalar exponential::H(const scalar p, const scalar expScalar, const scalar coeffSScalar) const
        {
            if (p < 0)
            {
                return 1.0 - coeffSScalar * pow(max(mag(p),VSMALL), expScalar);
            }
            else
            {
                return 1.0;
            }
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        exponential::exponential(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField)
            : saturationLaw(name, poroHydraulicProperties, pField),
              exponentialCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
              exp_(
                  IOobject(
                      "exp_swcc",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  mesh(),
                  dimensionedScalar(exponentialCoeffs_.lookup("exp"))),
              coeffS_(
                  IOobject(
                      "coeff_S",
                      mesh().time().timeName(),
                      db(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  mesh(),
                  dimensionedScalar(exponentialCoeffs_.lookup("coeffS"))),
              coeffk_(
                  IOobject(
                      "coeff_kr",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  mesh(),
                  dimensionedScalar(exponentialCoeffs_.lookup("coeffk"))),
              S_0(
                  IOobject(
                      "S_0",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  mesh(),
                  dimensionedScalar(exponentialCoeffs_.lookup("S_0")))

        {
            if (debug)
            {
                S_0.write();
                coeffk_.write();
                coeffS_.write();
                exp_.write();
            }
        }

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        scalar exponential::pStar(const label cellI) const
        {
            return GREAT;
        };

        scalar exponential::C(const scalar p, const label cellI)
        {
            if (p < 0)
            {
                scalar exp = exp_.primitiveField()[cellI];
                scalar coeffS = coeffS_.primitiveField()[cellI];
                return (exp * coeffS * pow(mag(p), exp)) / max(mag(p), SMALL);
            }
            else
            {
                return 0.0;
            }
        };

        scalar exponential::S(const scalar p, const label cellI)
        {
            scalar exp = exp_.primitiveField()[cellI];
            scalar coeffS = coeffS_.primitiveField()[cellI];
            return H(p,exp,coeffS) * S_0.internalField()[cellI];
        };

        scalar exponential::S
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            scalar exp = exp_.boundaryField()[patchI][faceI];
            scalar coeffS = coeffS_.boundaryField()[patchI][faceI];
            return H(p,exp,coeffS) * S_0.boundaryField()[patchI][faceI];
        };

        scalar exponential::kr(const scalar p, const label cellI)
        {
            scalar exp = exp_.primitiveField()[cellI];
            scalar coeffS = coeffS_.primitiveField()[cellI];
            scalar Hkr = H(p,exp,coeffS);
            scalar coeffk = coeffk_.primitiveField()[cellI];
            return (1 - coeffk * (1 - Hkr));
        };

        scalar exponential::kr
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            scalar exp = exp_.boundaryField()[patchI][faceI];
            scalar coeffS = coeffS_.boundaryField()[patchI][faceI];
            scalar Hkr = H(p,exp,coeffS);
            scalar coeffk = coeffk_.boundaryField()[patchI][faceI];
            return (1 - coeffk * (1 - Hkr));
        };



    } // End of namespace saturationLaws
} // End of namespace Foam

//*********************************************************** //
