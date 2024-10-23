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
    varSatPoroFluid

Description
    Richard's equation solver for flow through variably saturated porous media.

    This solver is based on p_rgh (pressure deviation from inital hydrostatic conditions)
    as the main variable.

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
    varSatPoroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroFluid_H
#define varSatPoroFluid_H

#include "variablySaturatedPoroFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace poroFluidModels
    {
        /*---------------------------------------------------------------------------*\
                        Class varSatPoroFluid Declaration
\*---------------------------------------------------------------------------*/

        class varSatPoroFluid
            : public variablySaturatedPoroFluid
        {
            // Private data

            //- Storage coefficient
            volScalarField Ss_;
            
            //- effective hydraulic conductivity devided by water specific weight and interpolated to surfaces
            autoPtr<surfaceScalarField> kEffbyGammafPtr_;

            bool steadyState_;

            // Private Member Functions

            //- Disallow default bitwise copy construct
            varSatPoroFluid(const varSatPoroFluid &);

            //- Disallow default bitwise assignment
            void operator=(const varSatPoroFluid &);

            void makeKEffbyGammaf();

            virtual void makePoroHydraulic();

        protected:

            surfaceScalarField &kEffbyGammaf();

        public:
            //- Runtime type information
            TypeName("varSatPoroFluid");
            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            varSatPoroFluid(
                Time &runTime,
                const word &region = dynamicFvMesh::defaultRegion,
                const bool sharedMesh = false);

            // Destructor

            virtual ~varSatPoroFluid()
            {
            }

            // Selector

            // Member Functions

            // Standard Access
            //- Renaming establisches p_rgh as pressure field
            virtual volScalarField &pField()
            {
                return p_rgh();
            }
                
            //- Renaming establisches p_rgh as pressure field
            virtual const volScalarField &pField() const
            {
                return p_rgh();
            }

            // Edit

            //- Evolve the fluid model
            virtual bool evolve();

            //- Return the desired new time-step size
            virtual scalar newDeltaT();

            virtual void writeFields(const Time &runTime);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace varSatPoroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
