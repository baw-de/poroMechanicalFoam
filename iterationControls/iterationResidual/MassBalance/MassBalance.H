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
    Foam::MassBalance

Description
    Using the total mass balance as convergence critera.

SeeAlso
    fvSolution

SourceFiles
    MassBalance.C

\*---------------------------------------------------------------------------*/

#ifndef delta_H
#define delta_H

#include "iterationResidual.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MassBalance Declaration
\*---------------------------------------------------------------------------*/

class MassBalance
:
    public iterationResidual
{
    // Private Member Functions
        bool writeField_;
        tmp<volScalarField> MassBalance_;
        HashTable<const poroFluidModel*> models

        void makeMassBalanceRef();

        //- No copy construct
        MassBalance(const MassBalance&) = delete;

public:
    //- runTime type information
        TypeName("MassBalance");

    // Constructors
        MassBalance
        (
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField
        );


    //- Destructor
    ~MassBalance();


    // Member Functions

        //- Execute specified function names
        virtual scalar calcResidual();

        virtual void reset();
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
