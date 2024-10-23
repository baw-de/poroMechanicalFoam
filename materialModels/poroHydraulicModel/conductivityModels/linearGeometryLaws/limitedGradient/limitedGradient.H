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
    Foam::conductivityModels::limitedGradient

Description
    Constitutive relations for soils: k(i), where k is greatly increased at icrit to
    limit the gradient to icrit.

Author
    Denis Maier, BAW. All rights reserved.
SourceFiles
    limitedGradient.C

\*---------------------------------------------------------------------------*/

#ifndef limitedGradient_H
#define limitedGradient_H

#include "conductivityModel.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace conductivityModels
    {

        /*---------------------------------------------------------------------------*\
                           Class limitedGradient Declaration
\*---------------------------------------------------------------------------*/

        class limitedGradient
            : public conductivityModel
        {
            // Private data
            word name_;
            dictionary& liquefyingZonePermeabilityProperties;
            dimensionedScalar kappa_;
            dimensionedScalar iCrit_;
            const volVectorField &i_;
            const dimensionedVector gamma_;
            volScalarField k0_;

            scalar kFunc(const scalar k0,const scalar iz) const;
            
        public:
            //- Runtime type information
            TypeName("limitedGradient");

            // Constructors

            //- Construct from components
            limitedGradient(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~limitedGradient()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual word name() { return "limitedGradient"; }

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
