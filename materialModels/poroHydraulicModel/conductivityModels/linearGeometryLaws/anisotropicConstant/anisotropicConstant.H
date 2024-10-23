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
    Foam::conductivityModels::anisotropicConstant

Description
    Expresses the cell anisotropic hydraulic conductivity tensor as scalar face orthogonal hydraulic conductivities
    on the cell faces.  

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    anisotropicConstant.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicConstant_H
#define anisotropicConstant_H

#include "conductivityModel.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace conductivityModels
    {

        /*---------------------------------------------------------------------------*\
                           Class anisotropicConstant Declaration
\*---------------------------------------------------------------------------*/

        class anisotropicConstant
            : public conductivityModel
        {
            // Private data
            word name_;
            dictionary& anisotropicConstantProperties;
            const dimensionedSymmTensor anisotropyFactors_;
            autoPtr<surfaceScalarField> kf_;       

        public:
            //- Runtime type information
            TypeName("anisotropicConstant");

            // Constructors

            //- Construct from components
            anisotropicConstant(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~anisotropicConstant()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary
            virtual word name() { return "anisotropicConstant"; }

            virtual tmp<surfaceScalarField> kf() const;

            virtual void updateK(){}

            virtual void write() const{};

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        };
    }; // End namespace conductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
