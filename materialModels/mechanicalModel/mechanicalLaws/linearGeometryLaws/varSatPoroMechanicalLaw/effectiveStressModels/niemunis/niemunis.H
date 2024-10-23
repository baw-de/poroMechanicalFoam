/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::niemunis

Description
    Turbulence model for niemunis incompressible flow.

SourceFiles
    niemunis.C

\*---------------------------------------------------------------------------*/

#ifndef niemunis_H
#define niemunis_H

#include "effectiveStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace effectiveStressModels
{
    /*---------------------------------------------------------------------------*\
                           Class niemunis Declaration
\*---------------------------------------------------------------------------*/

    class niemunis
        : public effectiveStressModel
    {

    public:
        //- Runtime type information
        TypeName("niemunis");

        // Constructors

        //- Construct from componentss
        niemunis(
            const dictionary &dict,
            const word effectiveStressModelName,
            const fvMesh &mesh);

        //- Destructor
        ~niemunis()
        {
        }

        // Member Functions

        //- Return 1 for therzaghi effective stress
        virtual tmp<scalarField> chi(const scalarField &n, const scalarField &S, const scalarField &p);

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace effectiveStressModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
