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
    Foam::richardsLinearization::Celia

Description
    This scheme implements modified Picard-method for the solution of Richard's equation.

    Reference:
    'A General Mass-Conservative Numerical-Solution for the Unsaturated Flow Equation'
    (Celia et al., 1990)

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
   Celia.C

\*---------------------------------------------------------------------------*/

#ifndef Celia_H
#define Celia_H

#include "richardsLinearization.H"
#include "varSatPoroHydraulicModel.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace richardsLinearizations
    {

        /*---------------------------------------------------------------------------*\
                           ClassCelia Declaration
\*---------------------------------------------------------------------------*/

        class Celia
            : public richardsLinearization
        {
        protected:
            // protected data
        public:
            //- runTime type information
            TypeName("Celia");

            // Constructors

            //- Construct from components
            Celia(
                const word &name,
                varSatPoroHydraulicModel &poroHydraulic,
                dictionary &poroFluidProperties,
                volScalarField& S);

            // Destructor

            ~Celia()
            {
            }

            // Member Functions

            virtual const word name() const { return name_; }

        virtual void initalize(volScalarField &totalP,volScalarField &pField);

        virtual bool checkConvergedAndUpdate(volScalarField &totalP, volScalarField &pField);

        virtual tmp<fvScalarMatrix> ddtS(const volScalarField &S, volScalarField &pField);

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace richardsLinearizations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
