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
    Foam::conductivityModels::gradientFunction

Description
    Constitutive relations for soils: k(i) 

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    gradientFunction.C

\*---------------------------------------------------------------------------*/

#ifndef gradientFunction_H
#define gradientFunction_H

#include "conductivityModel.H"
#include "volFields.H"
#include "Function1.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace conductivityModels
    {

        /*---------------------------------------------------------------------------*\
                           Class gradientFunction Declaration
\*---------------------------------------------------------------------------*/

        class gradientFunction
            : public conductivityModel
        {
            // Private data
            word name_;
            dictionary& gradientFunctionProperties;
            autoPtr<Function1<scalar>> kIFunc_;
            const volVectorField &i_;    
            const uniformDimensionedVectorField gamma_;     

        public:
            //- Runtime type information
            TypeName("gradientFunction");

            // Constructors

            //- Construct from components
            gradientFunction(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~gradientFunction()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual word name() { return "gradientFunction"; }

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
