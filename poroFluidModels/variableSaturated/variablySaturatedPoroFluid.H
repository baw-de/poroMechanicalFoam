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
    variablySaturatedPoroFluid

Description
    Base class for:
    Richard's equation solver for flow through variably saturated porous media.

    This solver can be based on p_rgh (pressure deviation from inital hydrostatic conditions)
    or pHead (relative/total pressure in m watercolumn) as the main variable.

    The equation considers the flow of pore fluid (water + gas bubbles). As that the
    Saturation S is considered the Saturation of pore fluid!
    This is an important distinction particularly when considering Saturation dependent
    Storage Ss, since here the Water Saturation Sw is used!
    The flux is converted from total volumetric flux to pore Fluid flux by the use
    of the **effective** hydraulic conductivity (parametric formuli).
    Also any source terms should also be considered to induce pore fluid, so that we
    need to take into account that not all the volume is filled with pore fluid.
    It thus needs to be scaled by saturation S. 

    Storage, hydraulic conductivity and Saturation are runtime-selectable functions
    handled by the 'varSatPoroHydraulicModel'-Object.

Author
    Denis Maier, BAW.  All rights reserved.

SourceFiles
    variablySaturatedPoroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef variablySaturatedPoroFluid_H
#define variablySaturatedPoroFluid_H

#include "poroFluidModel.H"
#include "varSatPoroHydraulicModel.H"
#include "richardsLinearization.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace poroFluidModels
    {
        /*---------------------------------------------------------------------------*\
                        Class variablySaturatedPoroFluid Declaration
\*---------------------------------------------------------------------------*/

        class variablySaturatedPoroFluid
            : public poroFluidModel
        {
            // Private data

            //- Disallow default bitwise copy construct
            variablySaturatedPoroFluid(const variablySaturatedPoroFluid &);

            //- Disallow default bitwise assignment
            void operator=(const variablySaturatedPoroFluid &);

        protected:
            // Protected Data

            //- Transport model
            mutable autoPtr<varSatPoroHydraulicModel> poroHydPtr_;

            //- Porohydraulic algorithm
            mutable autoPtr<richardsLinearization> richardsLinearizationPtr_;

            //- Saturation of pore fluid (Water + Bubbles)
            autoPtr<volScalarField> SPtr_;

            //- Mass Balance Field (demand driven)
            autoPtr<volScalarField> MassBalancePtr_;

            // Protected Member Functions
            virtual void makePoroHydraulic() = 0;

            void makeS();

            void makeRichardsLinearization();

            volScalarField& S()
            {
                if(!SPtr_.valid())
                {   
                    makeS();
                }
                return SPtr_.ref();
            }

            varSatPoroHydraulicModel& poroHydraulic()
            {
                if(!poroHydPtr_.valid())
                {   
                    makePoroHydraulic();
                }
                return poroHydPtr_.ref();
            }

            richardsLinearization& linearization()
            {
                if(!richardsLinearizationPtr_.valid())
                {   
                    makeRichardsLinearization();
                }
                return richardsLinearizationPtr_.ref();
            }


        public:
            //- Runtime type information
            TypeName("variablySaturatedPoroFluid");
            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            variablySaturatedPoroFluid(
                const word &type,
                Time &runTime,
                const word &fieldName,
                const word &region = dynamicFvMesh::defaultRegion,
                const bool sharedMesh = false);

            // Destructor

            virtual ~variablySaturatedPoroFluid()
            {
            }

            // Selector

            // Member Functions

            // Special Acess

            const volScalarField& S() const
            {
                return S();
            }

            virtual tmp<volScalarField> n0()
            {
                return poroHydraulic().n0();
            }

            const varSatPoroHydraulicModel& poroHydraulic() const
            {
                return poroHydraulic();
            }

            // Edit

            // Reminder: This still needs to be implemented by children
            //- Evolve the fluid model
            //virtual bool evolve() = 0;

            // Reminder: This still needs to be implemented by children
            //- Return the desired new time-step size
            //virtual scalar newDeltaT() = 0;

            void checkMassBalance();

            virtual scalar materialResidual()
            {
                return MassBalancePtr_.valid()
                           ? gSum(MassBalancePtr_())
                           : 0.0;
            }

            //- Return number of top-level iterations
            virtual label outerIters()
            {
                notImplemented("outerIters()");
                return iterCtrl().index(); // for now, to keep the compiler happy
            }

            // Reminder: This still needs to be implemented by children
            //virtual void writeFields(const Time &runTime) = 0;

            virtual void end();
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace variablySaturatedPoroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
