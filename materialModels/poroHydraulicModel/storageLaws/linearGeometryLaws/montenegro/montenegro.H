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
    Foam::storageLaws::montenegro

Description
    Extention of KPrime law with approximation of the influence of gas bubble compaction.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    montenegro.C

\*---------------------------------------------------------------------------*/

#ifndef montenegro_H
#define montenegro_H

#include "storageLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace storageLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class montenegro Declaration
\*---------------------------------------------------------------------------*/

        class montenegro
            : public storageLaw
        {
            // Private data
            dictionary& montenegroCoeffs_;
            const volScalarField S_pe_;
            const dimensionedScalar p_At_;
            const dimensionedScalar p_e_;

        public:
            //- Runtime type information
            TypeName("montenegro");

            // Constructors

            //- Construct from components
            montenegro(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~montenegro()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual scalar Ss
            (
                const scalar n,
                const scalar p,
                const label cellI
            );

            virtual bool write()
            {   
                return true;
            };

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace storageLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
