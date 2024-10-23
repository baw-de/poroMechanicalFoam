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
    Foam::iterationResidual

Description
    A class for calculating residualfields form solution fields.

SeeAlso
    fvSolution

SourceFiles
    iterationResidual.C

\*---------------------------------------------------------------------------*/

#ifndef iterationResidual_H
#define iterationResidual_H

#include "residualOperation.H"
#include "runTimeSelectionTables.H"
#include "Time.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class iterationResidual Declaration
\*---------------------------------------------------------------------------*/

class iterationResidual
{
    // Private Member Functions
        const Time& runTime_;
        const word name_;
        const scalar tolerance_;
        mutable autoPtr<residualOperation> operation_;

        //- No copy construct
        iterationResidual(const iterationResidual&) = delete;

        //- No copy assignment
        void operator=(const iterationResidual&) = delete;

protected:
        scalar residual_;
        const bool relative_;

        const Time& runTime()
        {
            return runTime_;
        }

public:

//- runTime type information
        TypeName("iterationResidual");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            iterationResidual,
            dictionary,
            (
                const Time& runTime,
                const word &name,
                const ITstream stream,
                const bool writeField),
            (runTime, name, stream, writeField));

        // Selectors

        //- Return a reference to the selected poroFluid model
        static autoPtr<iterationResidual> New(
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField);

        // Constructors

        //- Construct from components
        iterationResidual(
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField);


    //- Destructor
    virtual ~iterationResidual(){
        operation_.clear();
    };


    // Member Functions

    virtual scalar calcResidual() = 0;

    virtual word name() const
    {
        return name_;
    }

    scalar residual() const
    {
        return residual_;
    }

    scalar tolerance() const
    {
        return tolerance_;
    }

    scalar operation(const scalarField& x) const
    {
        return operation_->operation(x);
    }

    scalar operation(const List<scalar>& x) const
    {
        return operation_->operation(x);
    }

    word operationType() const
    {
        return operation_->name();
    }

    inline virtual void reset()
    {
        residual_ = GREAT;
    }

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
