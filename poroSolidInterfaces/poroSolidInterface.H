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
    poroSolidInterface

Description
    Virtual base class for poroFluid-solid interface coupling.
    We initialize solidModel object and a poroFluid object.
    For clarity of the code we will extract a reference
    the poroHydraulicModel from the poroFluid object.

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

SourceFiles
    poroSolidInterface.C

\*---------------------------------------------------------------------------*/

#ifndef poroSolidInterface_H
#define poroSolidInterface_H

#include "physicsModel.H"
#include "IOdictionary.H"
#include "dynamicFvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "poroFluidModel.H"
#include "solidModel.H"
#include "globalPolyPatch.H"
#include "meshToMesh.H"
#include "iterationControl.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                               Class poroSolidInterface Declaration
    \*---------------------------------------------------------------------------*/
    class poroSolidInterface
        : public physicsModel,
          public IOdictionary
    {
        // Private data

        //- reference to time
        Time& runTime_;

        //- FSI properties dictionary
        dictionary &interactionProperties_;

        //- Optional: is mesh shared between the two models?
        const Switch sharedMesh_;

        //- Optional: is porosity constant or updated?
        const Switch porosityConstant_;
        //- Optional: if porosity not constant: implicit or explicit treatment?
        const Switch porosityConstantExplicit_;

        //- Solid solver
        autoPtr<solidModel> solid_;

        //- poroFluid solver
        autoPtr<poroFluidModel> poroFluid_;

        //- Biot-Willis Coefficitent on the fluid mesh
        // It will be automatically grapped from the solid mesh
        // by the function checkMechanicalLawUpdateBiotCoeff
        autoPtr<volScalarField> b_;

        //- meshToMesh mapping
        autoPtr<meshToMesh> solidToPoroFluid_;

        //- List of residuals to be checked
        autoPtr<iterationControl> iterCtrl_;

        //- Internal work done by solid
        autoPtr<volScalarField> intWork_;

        // Private Member Functions

        void makeIterCtrl();

        void makePoroFluidCouplingSource();

        // Function to grap Biot-Willis Coefficitent from solid mesh
        bool checkMechanicalLawUpdateBiotCoeff
        (
            const mechanicalLaw& law,
            const label lawI,
            PtrList<volScalarField> &bs
        );

        void makeBiotCoeff();

        //- Disallow default bitwise copy construct
        poroSolidInterface(const poroSolidInterface &);

        //- Disallow default bitwise assignment
        void operator=(const poroSolidInterface &);

    protected:

        const Switch &porosityConstantExplicit() const
        {
            return porosityConstantExplicit_;
        }

    public:
        //- Runtime type information
        TypeName("poroSolidInterface");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            poroSolidInterface,
            dictionary,
            (
                Time & runTime,
                const word &region),
            (runTime, region));

        // Constructors

        //- Construct from components
        poroSolidInterface(
            const word &type,
            Time &runTime,
            const word &region = dynamicFvMesh::defaultRegion);

        // Selectors

        //- Select constructed from poroFluid and solid meshes
        static autoPtr<poroSolidInterface> New(
            Time &runTime,
            const word &region);

        // Destructor

        virtual ~poroSolidInterface();

    public:
        // Member Functions

        // Access
        //- Return const reference to the poroFluid main pressure field
        const volScalarField& pField() const
        {
            return poroFluid_().pField();
        }

        //- Return const reference to the poroFluid mesh
        const dynamicFvMesh &poroFluidMesh() const
        {
            return poroFluid_().mesh();
        }

        //- Return const reference to the solid mesh
        const dynamicFvMesh &solidMesh() const
        {
            return solid_().mesh();
        }

        //- Return time
        // to-do: consider storing a reference to time
        const Time &runTime() const
        {
            return runTime_;
        }

        //- Return sub-dictionary
        const dictionary &interactionProperties() const
        {
            return interactionProperties_;
        }

        //- Return poroFluid model
        poroFluidModel &poroFluidRef()
        {
            return poroFluid_();
        }

        //- Return poroFluid model
        const poroFluidModel &poroFluid()
        {
            return poroFluid_();
        }

        //- Return solid model
        solidModel &solidRef()
        {
            return solid_();
        }

        //- Return solid model
        const solidModel &solid()
        {
            return solid_();
        }

        //- Returns Biot-Willis coeff on poroFluid mesh
        const volScalarField b()
        {
            if(!b_.valid())
            {
                makeBiotCoeff();
            }
            return b_();
        }

        //- Return const access to the additional mesh correction switch
        const Switch &sharedMesh() const
        {
            return sharedMesh_;
        }

        //- Return const access to the additional porosity correction switch
        const Switch &porosityConstant() const
        {
            return porosityConstant_;
        }

        iterationControl& couplingControl()
        {
            if(!iterCtrl_.valid())
            {
                makeIterCtrl();
            }
            return iterCtrl_();
        }

        //- return reference to meshtoMesh interpolator
        const meshToMesh &solidToPoroFluid() const {
            return solidToPoroFluid_();
        }

        //- Calculate
        //- explicit coupling terms (This is virtual, meaning the function needs to be defined in
        //  children classes)
        virtual const tmp<volScalarField>  explicitCouplingDtoP() const = 0;
        //- explicit coupling terms (get multiplied with pField)
        virtual const tmp<volScalarField>  implicitCouplingDtoP() const = 0;

        // Edit

        //- Update the size of the time-step
        virtual void setDeltaT(Time &runTime);

        //- Evolve the interface
        virtual bool evolve() = 0;

        virtual void end();

        //- Initialize fields
        virtual void initializeFields();

        //- Update Total fields an potentially refine mesh
        void updateTotalFields();

        //- Move poroFluid mesh
        virtual void movePoroFluidMesh();

        //- Write fields
        virtual void writeFields(const Time &runTime);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //