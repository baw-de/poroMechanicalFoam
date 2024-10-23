/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the teL2 of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::L2

Description
    Sum of squares of resodual field. L2 Norm.

Author
    Denis Maier, BAW. All rights reserved.

SeeAlso
    fvSolution

SourceFiles
    L2.C

\*---------------------------------------------------------------------------*/

#ifndef L2_H
#define L2_H

#include "residualOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace residualOperations
    {

/*---------------------------------------------------------------------------*\
                         Class L2 Declaration
\*---------------------------------------------------------------------------*/

    class L2
    :
        public residualOperation
    {

            //- No copy construct
            L2(const L2&) = delete;

    public:
        //- runTime type information
            TypeName("L2");

        // Constructors
            L2
            (
                const word operation
            );


        //- Destructor
        ~L2();


        // Member Functions

            //- Execute specified function names
            virtual scalar operation(const scalarField& x) const;
            virtual scalar operation(const List<scalar>& x) const;

    };

} // End namespace residualOperations

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
