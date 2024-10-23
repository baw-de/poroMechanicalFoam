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

\*---------------------------------------------------------------------------*/

#include "terzaghi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace effectiveStressModels
    {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(terzaghi, 0);
        addToRunTimeSelectionTable(effectiveStressModel, terzaghi, dictionary);

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        terzaghi::terzaghi(
            const dictionary &dict,
            const word effectiveStressModelName,
            const fvMesh &mesh)
            : effectiveStressModel(dict, effectiveStressModelName, mesh)
        {
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        tmp<scalarField> terzaghi::chi(const scalarField &n, const scalarField &S, const scalarField &p)
        {
            tmp<scalarField> tchi
            (
                new scalarField(S.size(),1.0)
            );
            return tchi;
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace effectiveStressModels

} // End namespace Foam

// ************************************************************************* //
