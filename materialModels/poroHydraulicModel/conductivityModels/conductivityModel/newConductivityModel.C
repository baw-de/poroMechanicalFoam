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

#include "conductivityModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    autoPtr<conductivityModel> conductivityModel::New(
      const word &name,
      dictionary &poroHydraulicProperties,
      const volScalarField &pField)
    {

        Info << "Selecting hydraulic conductivity: "
             << name << endl;


#if (OPENFOAM >= 2112)
        auto* ctorPtr = dictionaryConstructorTable(name);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                name,
                "conductivityModel",
                name,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }
#else
        dictionaryConstructorTable::iterator cstrIter =
                dictionaryConstructorTablePtr_->find(name);

            if (cstrIter == dictionaryConstructorTablePtr_->end())
            {
                FatalErrorIn(
                    "conductivityModel::New(volScalarField&, "
                    "volScalarField&, "
                    "volScalarField&, "
                    "volScalarField&) ")
                    << "Unknown conductivityModel type "
                    << name << endl
                    << endl
                    << "Valid  conductivityModels are : " << endl
                    << dictionaryConstructorTablePtr_->toc()
                    << exit(FatalError);
            }

        auto* ctorPtr = cstrIter();
#endif
        return autoPtr<conductivityModel>(ctorPtr(name, poroHydraulicProperties, pField));
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
