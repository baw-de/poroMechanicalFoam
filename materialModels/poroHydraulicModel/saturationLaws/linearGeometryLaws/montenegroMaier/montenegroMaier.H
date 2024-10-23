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
    Foam::saturationLaws::montenegroMaier

Description
    soil constitutive relationship developed in-house at BAW. This relation incoporates
    soil storage from occluded gas bubbles into the SWCC. Therefor no storage model is needed. 

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    montenegroMaier.C

\*---------------------------------------------------------------------------*/

#ifndef montenegroMaier_H
#define montenegroMaier_H

#include "vanGenuchten.H"
#include "volFields.H"
#include "DimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class montenegroMaier Declaration
\*---------------------------------------------------------------------------*/

        class montenegroMaier
            : public vanGenuchten
        {
            // Private data

            dictionary& montenegroMaierCoeffs_;

            const volScalarField S_fgr0_;
            const dimensionedScalar p_At_;

            volScalarField Sf_;
            volScalarField Sfw_;
            volScalarField Ss_gr_;

            // Private Member Functions


            scalar Sfgr(const scalar p,const scalar Sfgr0) const;

            scalar Sfw(const scalar SfgrScalar) const;

            scalar Sf(const scalar &HSf, const scalar S0, const scalar Sr) const;

            scalar SFunc(const scalar SfwScalar, const scalar SfScalar) const;

            scalar Ss_gr(const scalar p, const label cellI) const;

        public:
            //- Runtime type information
            TypeName("montenegroMaier");

            // Constructors

            //- Construct from components
            montenegroMaier(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~montenegroMaier()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual scalar C(const scalar p, const label cellI);

            virtual scalar S(const scalar p, const label cellI);

            virtual scalar S
            (
                const scalar p,
                const label patchI,
                const label faceI
            );

            virtual void write() const{};

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace saturationLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
