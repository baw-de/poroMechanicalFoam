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
    poroFluidModel

Description
    Base class for solvers for flow through (variably saturated) porous media.

    Handles the initalization of basic field needed by every kind of flow solver
    and the accessability of these fields from coupling interfaces.

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
    poroFluidModel.C

\*---------------------------------------------------------------------------*/

#ifndef poroFluidModel_H
#define poroFluidModel_H

#include "physicsModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "dynamicFvMesh.H"
#include "fvc.H"
#include "fvOptions.H"
#include "iterationControl.H"
#include "uniformDimensionedFields.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        /*---------------------------------------------------------------------------*\
                                Class poroFluidModel Declaration
        \*---------------------------------------------------------------------------*/

        class poroFluidModel
            : public physicsModel, public dynamicFvMesh
        {
            //- Model name
            const word name_;

            //- Is the mesh a sharedMesh?
            const Switch sharedMesh_;

            //- reference to runTime
            const Time &runTime_;

            // Private data
            IOdictionary poroFluidProperties_;

            //- Derived additionalPhysicsModel type
            const word type_;

            //- name of the primary field of the additional physics
            const word fieldName_;

            //- gravity vector
            uniformDimensionedVectorField g_;

            //- pointer to the porosity field
            // It is either just a reference to the porosity in poroHydraulicModel
            // or, if porosity changes over time, its own field.
            autoPtr<volScalarField> nPtr_;

            // Initial porosity (Is only safed in case porosity gets updated)
            autoPtr<volScalarField> n0Ptr_;

            //- total pressure (actually relative to atmospheric pressure)
            volScalarField p_;

            //- main field time derivative
            volScalarField pDot_;

            //- pressure variation field header
            IOobject pRGHheader_;

            //- pressure variation from inital hydrostatic conditions
            volScalarField p_rgh_;

            //- main field gradient
            volVectorField hydraulicGradient_;

            // surface volumetric flux
            surfaceScalarField phi_;

            //- fvOption Sources, Constraints and Corrections
            fv::options& fvOptions_;

            //- name of the coeffs dictionary
            word coeffsDictName_;
            
            //- Control object for poroFluid loop
            autoPtr<iterationControl> iterCtrl_;

            //- maximal timestep when using adaptable timestep
            scalar maxDeltaT_;

        protected:
            //- Courant number when using adaptable timestep
            scalar CoNumber_;

            void makeN();

            //- This function will check that the D field was read from disk
            virtual void pRGHisRequired();

            //- Under-relaxation method
            //  It can be 'fixed', 'Aitken' or 'QuasiNewton'
            //const word relaxationMethod_;

        private:
            // Private Member Functions

            void makeIterCtrl();

            //- Disallow default bitwise copy construct
            poroFluidModel(const poroFluidModel &);

            //- Disallow default bitwise assignment
            void operator=(const poroFluidModel &);

        public:
            //- Runtime type information
            TypeName("poroFluidModel");

            // Declare run-time constructor selection table
            declareRunTimeSelectionTable(
                autoPtr,
                poroFluidModel,
                dictionary,
                (
                    Time & runTime,
                    const word &region,
                    const bool sharedMesh),
                (runTime, region, sharedMesh));

            // Constructors

            //- Construct from components
            poroFluidModel(
                const word &type,
                Time &runTime,
                const word &fieldName,
                const word &region = dynamicFvMesh::defaultRegion,
                const bool sharedMesh = false);
            // Destructor

            virtual ~poroFluidModel()
            {
            }

            // Selector

            static autoPtr<poroFluidModel> New(Time &runTime, const word &region, const bool sharedMesh = false);

            // Member Functions

            // Access

            //- return model name
            const word name() const
            {
                return name_;
            }

            bool sharedMesh() const
            {
                return sharedMesh_;
            }

            //- Return non-const mesh
            const dynamicFvMesh &mesh() const
            {
                return                 
                    sharedMesh_
                    ?runTime_.lookupObject<dynamicFvMesh>("solid")
                    :dynamic_cast<const dynamicFvMesh&>(*this);
            }

            //- Return time
            const Time &runTime() const
            {
                return runTime_;
            }

            dictionary &poroFluidDict()
            {
                return poroFluidProperties_.subDict(coeffsDictName_);
            }

            //- Return const reference to gravity vector
            const dimensionedVector& g() const
            {
                return g_;
            }

            virtual tmp<volScalarField> n0() = 0;

            //- Get the porosity field (This can either be the positity from poroHydraulicModel
            //  or a new porosity field generated to store changing porosity values)
            virtual const volScalarField& n(); 

            //- Update porosity
            virtual void update_porosity(const volScalarField Dn, const bool incremental);

            //- Return the main field (Might be any of p, pHead, p_rgh depending on solver)
            virtual volScalarField &pField() = 0;

            //- Return the main field (Might be any of p, pHead, p_rgh depending on solver)
            virtual const volScalarField &pField() const = 0;

            //- Return const reference to variable field
            volScalarField &p()
            {
                return p_;
            }

            //- Return const reference to variable field
            const volScalarField &p() const
            {
                return p_;
            }

            //- Return non-const reference to time derivative of variable field
            volScalarField &pDot()
            {
                return pDot_;
            }

            //- Return const reference to time derivative of variable field
            const volScalarField &pDot() const
            {
                return pDot_;
            }

            //- Return non-const reference to variable field
            volScalarField &p_rgh()
            {
                return p_rgh_;
            }

            //- Return const reference to variable field
            const volScalarField &p_rgh() const
            {
                return p_rgh_;
            }

            //- Return non-const reference to gradient of variable field
            volVectorField &hydraulicGradient()
            {
                return hydraulicGradient_;
            }

            //- Return const reference to gradient of variable field
            const volVectorField &hydraulicGradient() const
            {
                return hydraulicGradient_;
            }

            //- Return surface mass flux
            surfaceScalarField &phi()
            {
                return phi_;
            }

            //- Return const surface mass flux
            const surfaceScalarField &phi() const
            {
                return phi_;
            }

            fv::options& fvOptions()
            {
                return fvOptions_;
            }

            iterationControl& iterCtrl()
            {
                if(!iterCtrl_.valid())
                {
                    makeIterCtrl();
                }
                return iterCtrl_();
            }

            // Edit

            //- Evolve the fluid model
            virtual bool evolve() = 0;

            //- Return the desired new time-step size
            virtual scalar newDeltaT()
            {
                scalar maxDeltaTFact = CoNumber_*iterCtrl().nCycles()/(iterCtrl().index() + SMALL);
                scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
                return deltaTFact*runTime().deltaTValue();
            };

            virtual void setDeltaT(Time& runTime)
            {
                scalar newDT = min(newDeltaT(),maxDeltaT_);
                Info << "New automatic timestep: " << newDT << endl;
                runTime.setDeltaT
                (
                    newDT
                );
            }

            //- Return the mass balance residual
            virtual scalar materialResidual() = 0;

            //- Return number of top-level iterations
            virtual label outerIters() = 0;

            //- Write fields
            virtual void writeFields(const Time &runTime);

            //- End of the program
            virtual void end();

            //- Update or move poroFluid mesh
            virtual bool update();

            //- Write function must be defined for regIOobjects
            virtual bool writeData(Ostream &) const
            {
                return false;
            }

            //- calculate the L2-Norm of a scalarfield x
            scalar L2Norm(volScalarField x) const;
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
