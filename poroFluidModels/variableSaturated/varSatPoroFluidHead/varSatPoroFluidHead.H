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
    varSatPoroFluidHead

Description
    Richard's equation solver for flow through variably saturated porous media.

    This solver is based on pHead (relative/total pressure in m watercolumn)
    as the main variable.

    Storage, hydraulic conductivity and Saturation are runtime-selectable functions
    handled by the 'varSatPoroHydraulicModel'-Object.

Author
    Denis Maier, BAW.  All rights reserved.

SourceFiles
    varSatPoroFluidHead.C

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroFluidHead_H
#define varSatPoroFluidHead_H

#include "variablySaturatedPoroFluid.H"
#include "richardsLinearization.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
        namespace poroFluidModels
        {
        /*---------------------------------------------------------------------------*\
                        Class varSatPoroFluidHead Declaration
        \*---------------------------------------------------------------------------*/

        class varSatPoroFluidHead
            : public variablySaturatedPoroFluid
        {
            // Private data
            volScalarField pHead_;

            //- Storage coefficient
            volScalarField Ss_;

            //- effective hydraulic conductivity interpolated to surfaces
            autoPtr<surfaceScalarField> kEfffPtr_;

            bool steadyState_;

            // Private Member Functions

            //- Disallow default bitwise copy construct
            varSatPoroFluidHead(const varSatPoroFluidHead &);

            //- Disallow default bitwise assignment
            void operator=(const varSatPoroFluidHead &);

            void makeKEfff();

            virtual void makePoroHydraulic();

        protected:
            
            surfaceScalarField &kEfff();

        public:
            //- Runtime type information
            TypeName("varSatPoroFluidHead");

            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            varSatPoroFluidHead(
                    Time &runTime,
                    const word &region = dynamicFvMesh::defaultRegion,
                    const bool sharedMesh = false);

            // Destructor

            virtual ~varSatPoroFluidHead()
            {
            }

            // Selector

            // Member Functions

            // Access

                //- renameing of primary fields for convinience
            virtual volScalarField &pField()
            {
                return pHead_;
            }
                
                //- renameing of primary fields for convinience
            virtual const volScalarField &pField() const
            {
                return pHead_;
            }

            // Edit

            //- Evolve the fluid model
            virtual bool evolve();

            //- Return the desired new time-step size
            virtual scalar newDeltaT();

            virtual void writeFields(const Time &runTime);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        } // End namespace poroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
