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

\*---------------------------------------------------------------------------*/

#include "iterationResidual.H"
#include "deltaVf.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(iterationResidual, 0);
    defineRunTimeSelectionTable(iterationResidual, dictionary);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iterationResidual::iterationResidual
(
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField
)
:       runTime_(runTime),
        name_(name),
        // Tolerance can either be a number or the word "show"
        // in the last case the tolerance will be initiated as -1
        // which will always trigger convergence for this residual
        tolerance_(
            stream.last().isNumber()
            ? stream.last().number()
            : -1.0
        ),
        operation_(),
        residual_(GREAT),
        // bool if relative or total residual should be calculated
        relative_(stream.found(token(word("rel"))))
{
    relative_
         ? Info << "rel. "
         : Info << " ";
    Info  << name_
          << " Tol.: " <<  tolerance_ << endl;
    if (!stream.last().isNumber())
    {
        if(!(stream.last().wordToken() == "show"))
        {
            FatalErrorInFunction << "last entry for residual "
            << name_
            << "must be either a number or the word show"
            << ::Foam::abort(FatalError);
        }
    }
    operation_.reset(residualOperation::New(stream.first().wordToken()));

}

Foam::autoPtr<Foam::iterationResidual> Foam::iterationResidual::New(
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField)
{
    Info << "Initializing iteration control: " << nl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        return (autoPtr<Foam::iterationResidual>(new Foam::deltaVf(runTime, name, stream, writeField)));
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        return (autoPtr<Foam::iterationResidual>(new Foam::deltaVf(runTime, name, stream, writeField)));
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<Foam::iterationResidual>(ctorPtr(runTime, name, stream, writeField));

    // If not found try to make a iterational difference out of it

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Ostream operation  * * * * * * * * * * * * * * //

// ************************************************************************* //
