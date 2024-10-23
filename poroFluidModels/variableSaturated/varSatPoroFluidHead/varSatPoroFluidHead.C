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

#include "varSatPoroFluidHead.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
        namespace poroFluidModels
        {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(varSatPoroFluidHead, 0);
        addToRunTimeSelectionTable(poroFluidModel, varSatPoroFluidHead, dictionary);

        // * * * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

        void varSatPoroFluidHead::makeKEfff()
        {
            kEfffPtr_.reset(
                new surfaceScalarField(
                IOobject(
                      "makeKEfff",
                      runTime().timeName(),
                      pHead_.db(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  pHead_.mesh(),
                  dimensionedScalar("", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)
            ));
        }
    
        void varSatPoroFluidHead::makePoroHydraulic()
        {
            poroHydPtr_.reset(new varSatPoroHydraulicModel(pHead_, g()));
        }

        // * * * * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * * * //

        surfaceScalarField& varSatPoroFluidHead::kEfff()
        {
            if (!kEfffPtr_.valid())
            {
                makeKEfff();
            }

            return kEfffPtr_.ref();
        }



        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        varSatPoroFluidHead::varSatPoroFluidHead(
                Time &runTime,
                const word &region,
                const bool sharedmesh)
            : variablySaturatedPoroFluid(typeName, runTime, "pHead", region, sharedmesh),
              pHead_(
                IOobject(
                        "pHead",
                        runTime.timeName(),
                        // if mesh is shared between solid and addPhysics register 
                        // fields in this objectRegistery instead of (solid-) mesh
                        *this,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE),
                    mesh()),
              Ss_ //- Storage coefficient
              (
                (
                    IOobject
                    (
                        "Ss", 
                        runTime.timeName(),
                        pHead_.db(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    poroHydraulic().Ss(n(),pHead_)
                )
              ),
              kEfffPtr_(),
              steadyState_(word(mesh().ddtScheme("ddt(C,pHead)"))=="steadyState")
        {
            if 
            (
                !steadyState_ 
                && mesh().relaxField("kEfff")
                && mesh().fieldRelaxationFactor("kEfff")!=1.0
            )
            {
                FatalErrorIn("varSatPoroFluidHead::varSatPoroFluidHead") 
                    << " k should only be relaxed for steady-state calculations!!!"
                    << endl;
            }

            p() = pHead_ * poroHydraulic().magGamma();

            // Now that p is initalized we can initalize these object:
            SPtr_.reset(poroHydraulic().S(pHead_).ptr());

            // Initialize iteration control
            iterCtrl();
            // Check if iteration control checks mass balance,
            // if yes, we need to provide this field
            checkMassBalance();

            if (debug)
            {
                Info << mesh().sortedNames()
                     << endl;
            }

            Info << "Done generating variable saturated poroFluidModel (pHead)" << endl;
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        bool varSatPoroFluidHead::evolve()
        {

            // Initialize performance tracking
            // and silence automatic solver output

            SolverPerformance<scalar> solverPerfp;
            SolverPerformance<scalar>::debug=0;

            // Resetting all residuals
	        iterCtrl().reset();

            Info << "Evolving fluid model: " << this->type() << endl;
 
///////////- Save Fields from previous timeStep ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Update saturation
            S() = poroHydraulic().S(pHead_);
            // Update the effective hydraulic conductivity
            kEfff() = poroHydraulic().kEfff(pHead_);

            // Check if poroHydraulicModel changes storage
            // If yes, update storage
            if(poroHydraulic().updatesSs())
            {
                Ss_ = poroHydraulic().Ss(n(),pHead_);
            }

            ///////////////////- Update Flux ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Update hydraulic gradient for use in boundary conditions and models
            hydraulicGradient() = fvc::grad(pHead_) + fvc::grad(poroHydraulic().z());      
            // Update volumetric flux for use in boundary conditions
            phi() = -kEfff() * (mesh().magSf() * (fvc::snGrad(pHead_) + fvc::snGrad(poroHydraulic().z())));// Initialize flux phi() for BC that need the flux


///////////- Nonlinear Iterations (Conductivity+Boundary Condition, Picard/Fixed Point) ////////////////////////////////////////////////////////////////////////////////////////////////////
            do
            {
            	// Clear previous outer iterations convergence data
		        const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();

                // Store prev fields for relaxation
                pHead_.storePrevIter();

                // Initialize the iteration scheme
                linearization().initalize(pHead_,pHead_);
                //pHead_.correctBoundaryConditions(); 

//////////////////- Nonlinear Iterations (Saturation)  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                do
                {

                     if (debug)
		            {
		                Info << "Assembling pEqn"
		                    << endl;
		            
		                for (fv::option& source : static_cast<fv::optionList&>(fvOptions())){
		                    Info << "Source Terms " << source.name() << " is active? " << source.active() << endl; 
		                }
		                SolverPerformance<scalar>::debug=debug;
		            }

                // Assemble pressure equation
                fvScalarMatrix pEqn(                    
                    n() * linearization().ddtS(S(),pHead_)  //- Change in Saturation
                    + S() * fvm::ddt(Ss_, pHead_) //- Storage (e.g. compressibility of fluid)
                    - fvm::laplacian(kEfff(), pHead_)                       // pressure flux                                 
                    ==
                    fvc::laplacian(kEfff(), poroHydraulic().z())              // gravity flux
                    + fvOptions()(Ss_,pHead_)// optional: coupling and sources)
                );

                fvOptions().constrain(pEqn);   

                if (MassBalancePtr_.valid())
                {   
                    // We are taking the mass balance BEFORE solving, so we get the mass balance 
                    // from the inital pressure field. This means we are lagging behind 
                    // one iteration, but we include the nonlinearities form the coefficients
                    // that we otherwise wouldnt get.
                    MassBalancePtr_.ref().primitiveFieldRef() = pEqn.residual();
                }

                //Solve System
		        solverPerfp = Foam::solve(pEqn);
                
                fvOptions().correct(pHead_); 

                phi() = pEqn.flux();

                }while(!linearization().checkConvergedAndUpdate(pHead_,pHead_));

                // Relax pHead and kr (if steady state)
                pHead_.relax();                
                hydraulicGradient() = fvc::grad(pHead_) + fvc::grad(poroHydraulic().z());   

//////////////////- Output first and latest linear solver performance  ////////////////////////////////////////////////////////////////////////////////////////////////
                autoPtr<List<SolverPerformance<scalar>>> sp(new List<SolverPerformance<scalar>>(mesh().solverPerformanceDict().findEntry("pHead")->stream()));
                if(sp().size()==1)
                {
                    sp().last().print(Info.masterStream(mesh().comm()));
                }
                else
                {    
                    Info << "Initial solve:" << endl;
                    sp().first().print(Info.masterStream(mesh().comm()));
                    Info << "Final solve:" << endl;
                    sp().last().print(Info.masterStream(mesh().comm()));
                }

 ///////////////////- Update the coefficients //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Update saturation
                S() = poroHydraulic().S(pHead_);
                // Update the effective hydraulic conductivity
                
                // Check if poroHydraulicModel changes storage
                // If yes, update storage
                if(poroHydraulic().updatesSs())
                {
                    Ss_ = poroHydraulic().Ss(n(),pHead_);
                }

                if (steadyState_)
                {
                    kEfff().storePrevIter();
                    kEfff() = poroHydraulic().kEfff(pHead_);
                    kEfff().relax();
                }
                else
                {
                    kEfff() = poroHydraulic().kEfff(pHead_);

                }      

                pHead_.correctBoundaryConditions();

                }while (iterCtrl().loop());

                iterCtrl().write();
            
 ///////////////////- Update the coefficients after convergence //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            Info << "varSatPoroFluidHead evolved" << endl;
            p() = pHead_ * poroHydraulic().magGamma();
            p_rgh() = p() - poroHydraulic().p_Hyd();    
            //prgh_.correctBoundaryConditions();  
            pDot() = fvc::ddt(p());
            return 0;
        }

        scalar varSatPoroFluidHead::newDeltaT()
        {
            scalar maxCv = max(kEfff() / fvc::interpolate(Ss_+linearization().C(),"interpolate(Ss)")).value();
            Info << "max c_v: " << maxCv << endl;
            scalar minDeltaX = min(mesh().V() / fvc::surfaceSum(mag(mesh().Sf())));
            scalar wishedTimeStep = CoNumber_ * pow(minDeltaX, 2) / maxCv;
            scalar maxDeltaTFact = wishedTimeStep / runTime().deltaTValue();
            scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);
            return deltaTFact * runTime().deltaTValue();
        }

        void varSatPoroFluidHead::writeFields(const Time &runTime)
        {
            Info << "writing PhysicalModel Fields" << endl;

            volScalarField theta
            (
                "theta",
                S() * n()
            );
            theta.write();

            volScalarField kr
            (
                "kr",
                poroHydraulic().kr(p())
            );
            kr.write();

            volScalarField kEff
            (
                "kEff",
                kr*poroHydraulic().k()
            );
            kEff.write();

            volVectorField q(
                "q",
                fvc::reconstruct(phi())
            );
            q.write();

            poroFluidModel::writeFields(runTime);

        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace poroFluidModels
} // End namespace Foam

// ************************************************************************* //
