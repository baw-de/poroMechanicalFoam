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
    Foam::richardsLinearization::Casulli

Description
    This scheme implements Casulli-method for the solution of Richard's equation.

    Reference:
    Casulli & Zanolli, 2010:
    'A Nested Newton-Type Algorithm for Finite Volume Methods Solving
    Richards' Equation in Mixed Form'

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
   Casulli.C

\*---------------------------------------------------------------------------*/

#ifndef Casulli_H
#define Casulli_H

#include "richardsLinearization.H"
#include "varSatPoroHydraulicModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace richardsLinearizations
    {

        /*---------------------------------------------------------------------------*\
                           ClassCasulli Declaration
\*---------------------------------------------------------------------------*/

        class Casulli
            : public richardsLinearization
        {
            // Private data
            dictionary &casulliDict_;
            const scalar r_max_;
            const scalar nCasInt_;
            const scalar nCasExt_;

            volScalarField C_pStar_;
            volScalarField S_pStar_;
            volScalarField pStar_;

            autoPtr<volScalarField> P_;
            autoPtr<volScalarField> Q_;
            autoPtr<volScalarField> S1_;
            autoPtr<volScalarField> S2_;
            autoPtr<volScalarField> f1_;
            autoPtr<volScalarField> f2_;
            autoPtr<volScalarField> pPrevInt_;
            autoPtr<volScalarField> pPrevExt_;

            int extCorr_;
            int intCorr_;
            int totCorr_;

            label storedTimeIndex_;

        public:
            //- runTime type information
            TypeName("Casulli");

            // Constructors

            //- Construct from components
            Casulli(
                const word &name,
                varSatPoroHydraulicModel &poroHydraulic,
                dictionary &poroFluidProperties,
                volScalarField& S);

            // Destructor

            ~Casulli()
            {
            }

            // Member Functions

            virtual void initalize(volScalarField &totalP,volScalarField &pField);

            virtual bool checkConvergedAndUpdate(volScalarField &totalP, volScalarField &pField);

            virtual tmp<fvScalarMatrix> ddtS(const volScalarField &S, volScalarField &pField);

            void update_pStar(const volScalarField &new_pStar);

            void S1(const volScalarField &pField, const volScalarField &S);

            void P(const volScalarField &pField, const volScalarField &C);

            void S2(const volScalarField &pField, const volScalarField &S1Ref);

            void Q(const volScalarField &pField, const volScalarField &PRef);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace richardsLinearizations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
