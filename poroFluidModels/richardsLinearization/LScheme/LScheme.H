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
    Foam::richardsLinearization::LScheme

Description
    This scheme implements L-method for the solution of Richard's equation.
    
    Reference:
    'A study on iterative methods for solving Richards` equation' (List & Radu, 2016)

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
   LScheme.C

\*---------------------------------------------------------------------------*/

#ifndef LScheme_H
#define LScheme_H

#include "richardsLinearization.H"
#include "varSatPoroHydraulicModel.H"
#include "volFields.H"
#include "fvm.H"
#include "Celia.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace richardsLinearizations
    {

        /*---------------------------------------------------------------------------*\
                           ClassLScheme Declaration
\*---------------------------------------------------------------------------*/

        class LScheme
            : public Celia
        {
        public:
            //- runTime type information
            TypeName("LScheme");

            // Constructors

            //- Construct from components
            LScheme(
                const word &name,
                varSatPoroHydraulicModel &poroHydraulic,
                dictionary &poroFluidProperties,
                volScalarField& S);

            // Destructor

            ~LScheme()
            {
            }

            // Member Functions

            // Does nothing
            virtual void initalize(volScalarField &totalP,volScalarField &pField){}

            virtual tmp<fvScalarMatrix> ddtS(const volScalarField &S, volScalarField &pField);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace richardsLinearizations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
