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
    poroSolid

Description
    Solver for fully saturated poro-elasto-plasticity.

    We start the timeloop in which we update the coupling terms
    and first evolve the poroFluid solver, then evolve the 
    solid solver.

    The coupling from solid to fluid works by adding a 
    fvOption (poroSolidToFluidCouplingSource Class)
    to the fvOption list of the poroFluid object and 
    provide functions (explicit/implicitCouplingDtoP)
    that the fvOption accesses and adds it to the fluid 
    function.

    The coulpling from fluid to solid is handled by the
    wrapping class (poroMechanicalLaw/varSatPoroMechanicalLaw)
    around the solid mechanicalLaws. All it needs is to be able 
    to find a scalar field named p or p_rgh on the solid mesh/registry

Author
    Denis Maier, BAW. All rights reserved.
    Héctor Mentenegro, BAW. All rights reserved.

SourceFiles
    poroSolid.C

\*---------------------------------------------------------------------------*/

#ifndef poroSolid_H
#define poroSolid_H

#include "poroSolidInterface.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace poroSolidInteractions
    {

        /*---------------------------------------------------------------------------*\
                        Class poroSolid Declaration
\*---------------------------------------------------------------------------*/

        class poroSolid
            : public poroSolidInterface
        {
            // Private data

            //- Field for stabilization procedure for groundwater equation
            autoPtr<volScalarField> fixedStressStabil_;
            //- porosity change per unit time, which couples deformation into the flow equation
            autoPtr<volScalarField> nDot_;
            //- flow induced by solid acceleration, which couples deformation into the flow equation
            autoPtr<surfaceVectorField> q_relAcc_;
            //- momentum transfer from fluid to solid
            autoPtr<volScalarField> pRghSolidMesh_;
            //- momentum transfer from fluid to solid
            autoPtr<volScalarField> pSolidMesh_;

            // Private Member Functions

            //- Pressure Eqn coupling terms
            //- porosity change per unit time
            tmp<volScalarField> nDot(const volScalarField& b, const volVectorField& U);

            //- fixed stress stablization term
            tmp<volScalarField> fixedStressStabil(const volScalarField& b, const volScalarField& impK);

            //- fluxes arising from differencial acceleration (usually not significant)
            tmp<surfaceVectorField> q_relAcc(const surfaceScalarField& kf, const volVectorField& U);

            //- Disallow default bitwise copy construct
            poroSolid(const poroSolid &);

            //- Disallow default bitwise assignment
            void operator=(const poroSolid &);

        protected:
            // Protected member functions

        public:
            //- Runtime type information
            TypeName("poroSolid");

            // Constructors
            //- Construct from components
            poroSolid(
                Time &runTime,
                const word &region = "poroFluid");

            // Destructor

            virtual ~poroSolid()
            {
            }

            // Member Functions

            // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();
            
            // Access

            //- explicit coupling terms
            virtual const tmp<volScalarField>  explicitCouplingDtoP() const;
            //- explicit coupling terms (get multiplied with pField)
            virtual const tmp<volScalarField>  implicitCouplingDtoP() const; 


            virtual void writeFields(const Time &runTime);

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // namespace solidModels

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
