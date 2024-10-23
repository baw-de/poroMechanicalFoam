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

Class
    poroFluid

Description
    Solver for groundwater equation in pressure deviation from inital hydrostatic conditions (p_rgh) form.
    This is equivalent to groundwater potential form, however in pressure units
    instead of hight of watercolumn.

    Storage and hydraulic conductivity are runtime-selectable functions
    handled by the 'poroHydraulicModel'-Object.

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
    poroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef poroFluid_H
#define poroFluid_H

#include "poroFluidModel.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        namespace poroFluidModels
        {
            /*---------------------------------------------------------------------------*\
                                    Class poroFluid Declaration
            \*---------------------------------------------------------------------------*/

            class poroFluid
                : public poroFluidModel
            {
                
                //- Transport model
                mutable autoPtr<poroHydraulicModel> poroHydPtr_;

                //- Storage coefficient
                volScalarField Ss_;

                //- effective hydraulic conductivity devided by water specific weight and interpolated to surfaces
                surfaceScalarField kbyGammaf_;

                // Private Member Functions
                poroHydraulicModel& poroHydraulic();

                //- Disallow default bitwise copy construct
                poroFluid(const poroFluid &);

                //- Disallow default bitwise assignment
                void operator=(const poroFluid &);

            public:
                //- Runtime type information
                TypeName("poroFluid");
                // Declare run-time constructor selection table

                // Constructors

                //- Construct from components
                poroFluid(
                    Time &runTime,
                    const word &region = dynamicFvMesh::defaultRegion,
                    const bool sharedMesh = false);

                // Destructor

                virtual ~poroFluid()
                {
                }

                // Selector

                // Member Functions

                // Access
                
                //- renameing of primary fields for convinience
                virtual volScalarField& pField()
                {
                    return p_rgh();
                }
                
                //- renameing of primary fields for convinience
                virtual const volScalarField& pField() const
                {
                    return p_rgh();
                }

                virtual tmp<volScalarField> n0()
                {
                    return poroHydraulic().n0();
                }

                // Special Access
                const poroHydraulicModel& poroHydraulic() const
                {
                    return poroHydPtr_();
                }

                // Edit

                //- Evolve the fluid model
                virtual bool evolve();

                //- Return the desired new time-step size
                virtual scalar newDeltaT();

                //- since eqns are linear solverResidual is already mass balance, so no materialResidual is needed.
                virtual scalar materialResidual()
                {
                    return 0.0;
                }

                virtual label outerIters()
                {
                    return iterCtrl().index();
                }

                virtual void writeFields(const Time &runTime);

                virtual void end();
            };

            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        } // End namespace poroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
