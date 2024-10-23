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
    Foam::storageLaws::skemptonB

Description
    Implicitly describes soil and pore fluid compressibilities.
    It uses the Skempton parameter B and the Matrix compressibility to recalculate the 
    Storage.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    skemptonB.C

\*---------------------------------------------------------------------------*/

#ifndef skemptonB_H
#define skemptonB_H

#include "storageLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace storageLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class skemptonB Declaration
\*---------------------------------------------------------------------------*/

        class skemptonB
            : public storageLaw
        {
            // Private data
            dictionary& skemptonBCoeffs_;
            volScalarField BParam_;
            volScalarField Km_;

        public:
            //- Runtime type information
            TypeName("skemptonB");

            // Constructors

            //- Construct from components
            skemptonB(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~skemptonB()
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

            virtual bool updatesSs() const
            {
                return false;
            }

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace storageLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
