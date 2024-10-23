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
    Foam::saturationLaws::vanGenuchten

Description
    soil constitutive relationship as in van Genuchten (1980).

    Reference:
    van Genuchten, M.T. (1980). "A closed-form equation for predicting the hydraulic
    conductivity of unsaturated soils.". Soil Science Society of America Journal,
    44(5): 892-898

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    vanGenuchten.C

\*---------------------------------------------------------------------------*/

#ifndef vanGenuchten_H
#define vanGenuchten_H

#include "saturationLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class vanGenuchten Declaration
\*---------------------------------------------------------------------------*/

        class vanGenuchten
            : public saturationLaw
        {
            // Private data

            dictionary& vanGenuchtenCoeffs_;

        protected:
            // Protected Member Functions

            const volScalarField n_;
            volScalarField m_;
            const volScalarField S_r;
            const volScalarField S_0;
            const volScalarField alpha_;

            scalar SFunc
            (
                const scalar HS,
                const scalar S0,
                const scalar Sr
            );
            
            scalar krFunc
            (
                const scalar p,
                const scalar Hk,
                const scalar mk
            );

            scalar H
            (
                const scalar p,
                const scalar nH,
                const scalar mH,
                const scalar alphaH
            );


        public:
            //- Runtime type information
            TypeName("vanGenuchten");

            // Constructors

            //- Construct from components
            vanGenuchten(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~vanGenuchten()
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

            //virtual scalar dkbydp(const scalar p, const label cellI);

            virtual void write() const{};

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace saturationLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
