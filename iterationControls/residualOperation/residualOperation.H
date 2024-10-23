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
    Foam::residualOperation

Description
    A class for calculating residuals from residualfields.

Author
    Denis Maier, BAW. All rights reserved.

SeeAlso
    fvSolution

SourceFiles
    residualOperation.C

\*---------------------------------------------------------------------------*/

#ifndef residualOperation_H
#define residualOperation_H

#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class residualOperation Declaration
\*---------------------------------------------------------------------------*/

class residualOperation
{
    // Private Member Functions

        const word name_;
        //- No copy construct
        residualOperation(const residualOperation&) = delete;
        //- Disallow default bitwise assignment
        void operator=(const residualOperation &);
public:

//- runTime type information
        TypeName("residualOperation");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            residualOperation,
            dictionary,
            (
                const word operation
            ),
            (operation));

        // Selectors

        //- Return a reference to the selected poroFluid model
        static autoPtr<residualOperation> New(
            const word operation);

        // Constructors

        //- Construct from components
        residualOperation(
            const word operation);


    //- Destructor
    virtual ~residualOperation(){};


    // Member Functions
    word name() const {return name_;}
    virtual scalar operation(const scalarField& x) const = 0;
    virtual scalar operation(const List<scalar>& x) const = 0;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
