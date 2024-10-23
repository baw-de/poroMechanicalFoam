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
    Foam::conductivityModels::kozenyCarman

Description
    Constitutive relations for soils: k(n) 

    References:
    J. Kozeny, "Ueber kapillare Leitung des Wassers im Boden." Sitzungsber Akad. Wiss., Wien, 136(2a): 271-306, 1927.
    P.C. Carman, "Fluid flow through granular beds." Transactions, Institution of Chemical Engineers, London, 15: 150-166, 1937.
    P.C. Carman, "Flow of gases through porous media." Butterworths, London, 1956.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    kozenyCarman.C

\*---------------------------------------------------------------------------*/

#ifndef kozenyCarman_H
#define kozenyCarman_H

#include "conductivityModel.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace conductivityModels
    {

        /*---------------------------------------------------------------------------*\
                           Class kozenyCarman Declaration
\*---------------------------------------------------------------------------*/

        class kozenyCarman
            : public conductivityModel
        {
            // Private data
            word name_;
            dictionary& kozenyCarmanProperties;
            const dimensionedScalar D50_; 
            const dimensionedScalar mu_;  

            scalar kFunc(const scalar nScalar) const;

        public:
            //- Runtime type information
            TypeName("kozenyCarman");

            // Constructors

            //- Construct from components
            kozenyCarman(
                const word &name,
                dictionary &fluidProperties,
                const volScalarField &pField);

            // Destructor

            ~kozenyCarman()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual word name() { return "kozenyCarman"; }

            virtual scalar k(const label cellI);

            virtual scalar k
            (
                const label patchI,
                const label faceI
            );

            virtual void write() const{};

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        };
    }; // End namespace conductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
