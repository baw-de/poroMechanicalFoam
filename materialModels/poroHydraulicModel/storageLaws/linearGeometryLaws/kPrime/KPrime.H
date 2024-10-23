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
    Foam::storageLaws::KPrime

Description
    Law, where storage is dependend of Saturation and (optionally) pressure.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    KPrime.C

\*---------------------------------------------------------------------------*/

#ifndef KPrime_H
#define KPrime_H

#include "storageLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace storageLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class KPrime Declaration
\*---------------------------------------------------------------------------*/

        class KPrime
            : public storageLaw
        {
            // Private data
            dictionary& KPrimeCoeffs_;
            // Water-Saturation at 0 Pressure
            // (not to be confused with fluid-saturation S_0)
            const volScalarField S_p0_;
            const dimensionedScalar Kw_;
            const dimensionedScalar p_At_;
            const Switch pDep_;

        public:
            //- Runtime type information
            TypeName("KPrime");

            // Constructors

            //- Construct from components
            KPrime(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~KPrime()
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
                if(pDep_ || writeSs_)
                {
                    writeSs_ = false;
                    return true;
                }
                return false;
            };

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace storageLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
