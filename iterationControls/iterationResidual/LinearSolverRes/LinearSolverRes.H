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
    Foam::LinearSolverRes

Description
    Using the linear solver inital residuals as convergence critera.

Author
    Denis Maier, BAW. All rights reserved.


SeeAlso
    fvSolution

SourceFiles
    LinearSolverRes.C

\*---------------------------------------------------------------------------*/

#ifndef delta_H
#define delta_H

#include "iterationResidual.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LinearSolverRes Declaration
\*---------------------------------------------------------------------------*/

class LinearSolverRes
:
    public iterationResidual
{
    // Private Member Functions
        //- No copy construct
        LinearSolverRes(const LinearSolverRes&) = delete;

public:
    //- runTime type information
        TypeName("LinearSolver");

    // Constructors
        LinearSolverRes
        (
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField
        );


    //- Destructor
    ~LinearSolverRes();


    // Member Functions

        //- Execute specified function names
        virtual scalar calcResidual();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
