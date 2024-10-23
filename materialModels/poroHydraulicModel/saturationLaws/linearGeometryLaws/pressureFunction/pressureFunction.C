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

#include "pressureFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {
        defineTypeNameAndDebug(pressureFunction, 0);

        addToRunTimeSelectionTable(
            saturationLaw,
            pressureFunction,
            dictionary);

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        pressureFunction::pressureFunction(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField)
            : saturationLaw(name, poroHydraulicProperties, pField),
              pressureFunctionCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
              pStar_(pressureFunctionCoeffs_.lookupOrAddDefault<dimensionedScalar>("CFuncMaximaPoint",dimensionedScalar("pStar",pField.dimensions(),-GREAT))),
              SFunc_(Function1<scalar>::New("SOfP",pressureFunctionCoeffs_.subDict("S"))),
              CFunc_(Function1<scalar>::NewIfPresent("CofP",pressureFunctionCoeffs_.subDict("C"))),
              krSFunc_(Function1<scalar>::NewIfPresent("krOfS",pressureFunctionCoeffs_.subDict("kr"))),
              krPFunc_(Function1<scalar>::NewIfPresent("krOfP",pressureFunctionCoeffs_.subDict("kr")))
        {
            if(krSFunc_.valid() && krPFunc_.valid())
            {
                FatalErrorInFunction() << "relative hydr. conductivity can only be a function of pressure OR saturation"
                                    << endl;
            }
            else if(!krSFunc_.valid() && !krPFunc_.valid())
            {
                FatalErrorInFunction() << "relative hydr. conductivity needs either pressure or saturation function"
                                    << endl;
            }
        }

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        scalar pressureFunction::pStar(const label cellI) const
        {
            return pStar_.value();
        }

        scalar pressureFunction::C(const scalar p, const label cellI)
        {
            if(!CFunc_.valid())
            {
                return (SFunc_().value(p+SMALL)-SFunc_().value(p-SMALL))/(2*SMALL);
            }
            else
            {
                return CFunc_().value(p);
            }
        }

        scalar pressureFunction::S(const scalar p, const label cellI)
        {
            return SFunc_().value(p);
        }

        scalar pressureFunction::S
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            return S(p,patchI);
        }

        scalar pressureFunction::kr(const scalar p, const label cellI)
        {
            if(!krSFunc_.valid())
            {
                return krPFunc_().value(p);
            }
            else
            {
                return krSFunc_().value(SFunc_().value(p));
            }
        }

        scalar pressureFunction::kr
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            return kr(p,patchI);
        }

    } // End of namespace saturationLaws
} // End of namespace Foam

//*********************************************************** //
