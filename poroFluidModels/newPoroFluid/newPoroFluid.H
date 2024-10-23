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
    newPoroFluid

Description
    Richard's equation solver for flow through variably saturated porous media

Author
    Denis Maier, BAW.  All rights reserved.

SourceFiles
    newPoroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef newPoroFluid_H
#define newPoroFluid_H

#include "poroFluidModel.H"
#include "poroHydraulicModel.H"
#include "richardsLinearization.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
        namespace poroFluidModels
        {
        /*---------------------------------------------------------------------------*\
                        Class newPoroFluid Declaration
\*---------------------------------------------------------------------------*/

        class newPoroFluid
            : public poroFluidModel
        {
            // Private data

            //- Transport model
            mutable autoPtr<poroHydraulicModel> poroHydPtr_;

            //- Porohydraulic algorithm
            mutable autoPtr<richardsLinearization> richardsLinearizationPtr_;

            surfaceScalarField kEffbyGammaf_;

            autoPtr<volScalarField> MassBalance_;

            bool steadyState_;

            // Private Member Functions
            void updateSaturation(const volScalarField& p);

            void updateFlux();

            //- return total pressure
            tmp<volScalarField> p();

            //- Disallow default bitwise copy construct
            newPoroFluid(const newPoroFluid &);

            //- Disallow default bitwise assignment
            void operator=(const newPoroFluid &);

        public:
            //- Runtime type information
            TypeName("newPoroFluid");
            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            newPoroFluid(
                    Time &runTime,
                    const word &region = dynamicFvMesh::defaultRegion,
                    const bool sharedMesh = false);


            // Destructor

            virtual ~newPoroFluid()
            {
            }

            // Selector

            // Member Functions

            // Access
                //- renameing of primary fields for convinience
                volScalarField &p_rgh()
                {
                    return vf();
                }

                volScalarField &p_rghDot()
                {
                    return vfDot();
                }
                //- renameing of primary fields for convinience
                const volScalarField &p_rgh() const
                {
                    return vf();
                }

                const volScalarField &p_rghDot() const
                {
                    return vfDot();
                }

                virtual const tmp<surfaceScalarField> k()
                {
                    tmp<surfaceScalarField> tk(
                    new surfaceScalarField(kEffbyGammaf_*poroHydraulic().magGamma()));
                    return tk;
                }

            // Edit

            //- Evolve the fluid model
            virtual bool evolve();

            //- Return the desired new time-step size
            virtual scalar newDeltaT()
            {
                return VGREAT;
            };

            virtual scalar materialResidual()
            {               
                return 
                MassBalance_.valid()
                ?gSum(MassBalance_())
                :0.0;
            }

            //- Return number of top-level iterations
            virtual label outerIters()
            {
                notImplemented("outerIters()");
                return iterCtrl().index(); //for now, to keep the compiler happy
            }

            virtual void writeFields(const Time &runTime);

            poroHydraulicModel &poroHydraulic() const;

            richardsLinearization &linearization() const;

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace poroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
