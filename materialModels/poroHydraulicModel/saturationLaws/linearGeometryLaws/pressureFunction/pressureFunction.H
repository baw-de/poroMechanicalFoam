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

Class
    Foam::saturationLaws::pressureFunction

Description
    Soil consistutive relations, where S(p) is a runtime defineable expression.
    Optionally, a relation C(p) for the exact jacobian can be difined as well, otherwise
    it will be calulated though central differencing.
    For relative hydraulic conductivity, kr can be a runtime definable expression of 
    saturation OR pressure.

Author
    Denis Maier, BAW. All rights reserved.


SourceFiles
    pressureFunction.C

\*---------------------------------------------------------------------------*/

#ifndef pressureFunction_H
#define pressureFunction_H

#include "saturationLaw.H"
#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class pressureFunction Declaration
\*---------------------------------------------------------------------------*/

        class pressureFunction
            : public saturationLaw
        {
            // Private data

            dictionary& pressureFunctionCoeffs_;
            dimensionedScalar pStar_;
            autoPtr<Function1<scalar>> SFunc_;
            autoPtr<Function1<scalar>> CFunc_;
            autoPtr<Function1<scalar>> krSFunc_;
            autoPtr<Function1<scalar>> krPFunc_;

            // Private Member Functions

        public:
            //- Runtime type information
            TypeName("pressureFunction");

            // Constructors

            //- Construct from components
            pressureFunction(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~pressureFunction()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual scalar pStar(const label cellI) const;

            virtual scalar C(const scalar p, const label cellI);

            virtual scalar S(const scalar p, const label cellI);

            virtual scalar S
            (
                const scalar p,
                const label patchI,
                const label faceI
            );

            virtual scalar kr(const scalar p, const label cellI);

            virtual scalar kr
            (
                const scalar p,
                const label patchI,
                const label faceI
            );

            virtual void write() const{};

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace saturationLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
