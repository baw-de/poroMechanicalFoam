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
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::max

Description
    Maximum of residal field. L_Inf norm.

Author
    Denis Maier, BAW. All rights reserved.

SeeAlso
    fvSolution

SourceFiles
    max.C

\*---------------------------------------------------------------------------*/

#ifndef maxRes_H
#define maxRes_H

#include "residualOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace residualOperations
    {

/*---------------------------------------------------------------------------*\
                         Class max Declaration
\*---------------------------------------------------------------------------*/

    class max
    :
        public residualOperation
    {

            //- No copy construct
            max(const max&) = delete;

    public:
        //- runTime type information
            TypeName("max");

        // Constructors
            max
            (
                const word operation
            );


        //- Destructor
        ~max();


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
