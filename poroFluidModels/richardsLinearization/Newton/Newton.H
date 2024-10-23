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
    Foam::richardsLinearization::Newton

Description
    This scheme implements Newton's method for the solution of Richard's equation.
    (Not compiled at the moment: bad performance)
    
    Reference:
    'A study on iterative methods for solving Richards` equation' (List & Radu, 2016)

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
   Newton.C

\*---------------------------------------------------------------------------*/

#ifndef Newton_H
#define Newton_H

#include "richardsLinearization.H"
#include "poroHydraulicModel.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace richardsLinearizations
    {

        /*---------------------------------------------------------------------------*\
                           ClassNewton Declaration
\*---------------------------------------------------------------------------*/

        class Newton
            : public richardsLinearization
        {
        protected:
            // protected data
        public:
            //- runTime type information
            TypeName("Newton");

            // Constructors

            //- Construct from components
            Newton(
                const word &name,
                poroHydraulicModel &poroHydraulic,
                dictionary &poroFluidProperties);

            // Destructor

            ~Newton()
            {
            }

            // Member Functions

            virtual const word name() const { return name_; }

        virtual void initalize(volScalarField &totalP,volScalarField &pField);

        virtual bool checkConvergedAndUpdate(volScalarField &totalP, volScalarField &pField);

        virtual tmp<fvScalarMatrix> ddtS(const volScalarField &S, volScalarField &pField);

        virtual tmp<fvScalarMatrix> ddpk(const surfaceScalarField &kField, volScalarField &pField);

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace richardsLinearizations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
