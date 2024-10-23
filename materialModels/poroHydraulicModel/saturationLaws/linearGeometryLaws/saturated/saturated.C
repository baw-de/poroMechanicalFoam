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

\*---------------------------------------------------------------------------*/

#include "saturated.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {
        defineTypeNameAndDebug(saturated, 0);

        addToRunTimeSelectionTable(
            saturationLaw,
            saturated,
            dictionary);

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        saturated::saturated(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField)
            : saturationLaw(name, poroHydraulicProperties, pField)
        {}

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        scalar saturated::pStar(const label cellI) const
        {
            return -GREAT;
        }

        scalar saturated::C(const scalar p, const label cellI)
        {
            return 0.0;
        }

        scalar saturated::S(const scalar p, const label cellI)
        {
            return 1.0;
        }

        scalar saturated::S
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            return 1.0;
        }

        scalar saturated::kr(const scalar p, const label cellI)
        {
            return 1.0;
        }

        scalar saturated::kr
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            return 1.0;
        }

    } // End of namespace saturationLaws
} // End of namespace Foam

//*********************************************************** //
