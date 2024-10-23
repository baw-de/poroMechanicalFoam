/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "variablySaturatedPoroFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace poroFluidModels
    {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(variablySaturatedPoroFluid, 0);

        // * * * * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * * * //

        void variablySaturatedPoroFluid::makeS()
        {
            SPtr_.reset
            (
                new volScalarField
                (
                        IOobject
                        (
                            "S", 
                            pField().time().timeName(),
                            pField().db(),
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        pField().mesh(),
                        dimensionedScalar(dimless,0.0),
                        "zeroGradient"
                )
            );
        }

        void variablySaturatedPoroFluid::makeRichardsLinearization()
        {   
            // Initialize the linearization method
            richardsLinearizationPtr_ =
                richardsLinearization::New
                (
                    poroFluidDict().get<word>("solutionAlgorithm"),
                    poroHydraulic(),
                    poroFluidDict(),
                    S()
                );
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        variablySaturatedPoroFluid::variablySaturatedPoroFluid
        (
            const word& type,
            Time& runTime,
            const word& fieldName,
            const word& region,
            const bool sharedMesh
        )
            : poroFluidModel(type, runTime, fieldName, region, sharedMesh),
              poroHydPtr_(), // Ptr to Unified saturation and storage models 
              richardsLinearizationPtr_(), // How to linearize Richard's equation
              SPtr_(), //We need to update the total pore pressure before initilizing Saturation
              MassBalancePtr_()
        {}

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        void variablySaturatedPoroFluid::checkMassBalance()
        {
            forAll(iterCtrl().residuals(), iRes)
            {
                if (iterCtrl().residuals()[iRes].name() == "MassBalance")
                {
                    MassBalancePtr_.reset(
                        new volScalarField(
                            IOobject(
                                "MassBalanceResidual",
                                pField().time().timeName(),
                                pField().db(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE),
                            pField().mesh(),
                            dimensionedScalar("", dimVolume, 0.0)));
                }
            }
        }

        void variablySaturatedPoroFluid::end()
        {
            poroHydraulic().write();
            poroFluidModel::end();
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace variablySaturatedPoroFluidModels
} // End namespace Foam

// ************************************************************************* //
