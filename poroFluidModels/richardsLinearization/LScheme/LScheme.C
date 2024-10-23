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

#include "LScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    namespace richardsLinearizations
    {
        defineTypeNameAndDebug(LScheme, 0);

        addToRunTimeSelectionTable(
            richardsLinearization,
            LScheme,
            dictionary);

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        LScheme::LScheme(
            const word &name,
            varSatPoroHydraulicModel &poroHydraulic,
            dictionary &poroFluidProperties,
            volScalarField& S)
            : Celia(name, poroHydraulic, poroFluidProperties,S)
        {
            C().rename("L");
            C() = poroHydraulic.C(poroHydraulic.pStar());
        }

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        tmp<fvScalarMatrix> LScheme::ddtS(const volScalarField &S, volScalarField &pField)
        {
            tmp<fvScalarMatrix> ddtSTMP( new fvScalarMatrix(
                pField,
                C().dimensions() * pField.dimensions() * dimVol / dimTime));

            scalar rDeltaT = 1.0 / mesh().time().deltaT().value();

            fvScalarMatrix &ddtS = ddtSTMP.ref();
            
            ddtS.diag() = rDeltaT * C().internalField() * mesh().V();
            ddtS.source() = rDeltaT * (S.oldTime() - S + C() * pField.prevIter()) * mesh().V();
            return ddtSTMP;
        }

    } // End of namespace richardsLinearizations
} // End of namespace Foam

//*********************************************************** //
