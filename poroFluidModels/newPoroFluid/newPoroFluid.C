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

#include "newPoroFluid.H"
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

        defineTypeNameAndDebug(newPoroFluid, 0);
        addToRunTimeSelectionTable(poroFluidModel, newPoroFluid, dictionary);

        // * * * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

        tmp<volScalarField> newPoroFluid::p(){
            return p_rgh() + poroHydraulic().p_Hyd();
        }

        void newPoroFluid::updateSaturation(const volScalarField& p) 
        {
            poroHydraulic().S(p);
            poroHydraulic().kr(p);
            if(steadyState_)
            {
                    const_cast<volScalarField&>(poroHydraulic().kr()).relax();
            }
            kEffbyGammaf_ = //fvc::interpolate(poroHydraulic().kr()) * poroHydraulic().k()/poroHydraulic().magGamma();
                            poroHydraulic().keff()/poroHydraulic().magGamma();
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        newPoroFluid::newPoroFluid(
                Time &runTime,
                const word &region,
                const bool sharedmesh)
            : poroFluidModel(typeName, runTime, "p_rgh", region, sharedmesh),
              poroHydPtr_(),    // Ptr to Unified saturation and storage models
              richardsLinearizationPtr_(), // Ptr to solution algorithm for groundwater eqn
              kEffbyGammaf_(                                           // effective hydraulic conductivity (saturated * unsaturated)/gamma_w
                  IOobject(
                      "kEffbyGammaf",
                      runTime.timeName(),
                      p_rgh().db(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                    p_rgh().mesh(),
                  dimensionedScalar("",dimensionSet(-1,3,1,0,0,0,0),0.0)),
              MassBalance_(),
              steadyState_(word(mesh().ddtScheme("ddt(C,p_rgh)"))=="steadyState")
        {
            //- Check if we have a steady state case (in which case underrelaxation for kr can be activated, otherwise not)
            if ( !steadyState_
                && mesh().relaxField("kr") && mesh().fieldRelaxationFactor("kr")!=1.0)
            {
                FatalErrorIn("newPoroFluidHead::newPoroFluidHead") << " kr should only be relaxed for steady-state calculations!!!"
                                                                   << endl;
            }  

            // Build poroHydraulicModel
            poroHydPtr_.reset(poroHydraulicModel::New(p_rgh(),g()));
            updateSaturation(p());
            poroHydraulic().Ss(p());
            richardsLinearizationPtr_ =
                  richardsLinearization::New(
                      word(poroFluidDict().lookup("solutionAlgorithm")),
                      poroHydraulic(),
                      poroFluidDict());

            iterCtrl();
            forAll(iterCtrl().residuals(), iRes)
            {
                if(iterCtrl().residuals()[iRes].name()=="MassBalance")
                {
                    MassBalance_.reset(
                        new volScalarField(
                        IOobject(
                        "MassBalanceResidual",
                        runTime.timeName(),
                        p_rgh().db(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                        p_rgh().mesh(),
                    dimensionedScalar("",dimVolume,0.0)));
                }
            }

            if (debug)
            {
                Info << p_rgh().db().sortedNames()
                     << endl;
            }

            Info << "Done generating extended poroFluidModel" <<endl;
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        bool newPoroFluid::evolve()
        {

            // Initialize performance tracking
            // and silence automatic solver output

            SolverPerformance<scalar> solverPerfp;
            SolverPerformance<scalar>::debug=0;


            Info << "Evolving fluid model: " << this->type() << endl;

            ///////////- Save Fields from previous timeStep ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            volScalarField p_ = p();
            updateSaturation(p_);
            poroHydraulic().Ss(p_);

            linearization(); // Initialize pEqn solver (linearization() --> makerichardsLinearization())
            
            ///////////////////- Update Flux ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
            hydraulicGradient() = fvc::grad(p_rgh())/poroHydraulic().magGamma();
            poroHydraulic().updateFlux(phi(), hydraulicGradient());
            volScalarField divQ("divQ",fvc::div(phi()));
            ///////////- Nonlinear Iterations (Flux/Conductivity+Boundary Condition) ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            iterCtrl().reset();
            do 
            {
            	// Clear previous outer iterations convergence data
		        const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();

                // Store prev fields for relaxation
                poroHydraulic().kr().storePrevIter();
                p_rgh().storePrevIter();
                poroHydraulic().S().storePrevIter();
                phi().storePrevIter();

                ///////////////////- Initialize/Update the coefficients //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                linearization().initalize(p_,p_rgh());
                p_rgh() = p_ - poroHydraulic().p_Hyd(); // In case Casulli's Scheme is used, p_ is capped initially, this needs to be transfered to p_rgh
                //p_rgh().correctBoundaryConditions();
                volScalarField explicitLaplacian = fvc::laplacian(kEffbyGammaf_, p_rgh());

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
		            
		            fvScalarMatrix pEqn(
		            		poroHydraulic().n() * linearization().ddtS(poroHydraulic().S(),p_rgh()) 	    //- Change in Saturation
		            		+ poroHydraulic().n() * fvm::ddt(poroHydraulic().Ss(), p_rgh())     	        //- Storage (e.g. compressibility of fluid)
                        	==
                            fvm::laplacian(kEffbyGammaf_, p_rgh()) 					                        //- Implicit darcian fluxes estimate
                            - explicitLaplacian 					                                        //- Explicit darcian fluxes estimate (must be calculated after saturation loop to ensure consistency of fluxes)
                            - divQ //fvc::div(phi())                                                               //- Fluxes
		            		+ fvOptions()(poroHydraulic().Ss(),p_rgh())					                    //- Source terms (e.g. volumetric deformation, wells etc.)
		            );

		            fvOptions().constrain(pEqn);

                    //Solve System
		            solverPerfp = pEqn.solve();  
		            	
		            fvOptions().correct(p_rgh());
		            p_ = p(); // Update the total pressure with new excess pressure
		            
		                if (debug)
		                {
		                    Info << "Solved pEqn" << endl;
		                }
		                
		        }while(!linearization().checkConvergedAndUpdate(p_,p_rgh())); //Check Convergence and Update Internal Veriables
                //////////////////- Nonlinear Iterations (Saturation) converged  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Relax p_rgh
                p_rgh().relax();
                p_ = p();

                //////////////////- Output first and latest linear solver performance  ////////////////////////////////////////////////////////////////////////////////////////////////
                autoPtr<List<SolverPerformance<scalar>>> sp(new List<SolverPerformance<scalar>>(mesh().solverPerformanceDict().findEntry("p_rgh")->stream()));
                if(sp().size()==1)
                {
                    sp().last().print(Info.masterStream(mesh().comm()));
                }
                else
                {   
                    Info << "Initial linear solver run:" << endl;
                    sp().first().print(Info.masterStream(mesh().comm()));
                    Info << "Final linear solver run:" << endl;
                    sp().last().print(Info.masterStream(mesh().comm()));
                }

                ///////////////////- Update the Coeffs and secondary fields ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
                updateSaturation(p_);
                poroHydraulic().Ss(p_);

                ///////////////////- Update Flux ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
                hydraulicGradient() = fvc::grad(p_rgh())/poroHydraulic().magGamma();
                poroHydraulic().updateFlux(phi(), hydraulicGradient());
                phi().relax();
                divQ = fvc::div(phi());
                p_rgh().correctBoundaryConditions();

                ///////////////////- Calculate Mass Balance Errors ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(MassBalance_.valid())
                {
                    MassBalance_.ref().primitiveFieldRef() = (poroHydraulic().n()*(fvc::ddt(poroHydraulic().S())+fvc::ddt(poroHydraulic().Ss(),p_rgh())) * mesh().V()
                                - fvc::div(phi())().primitiveField() * mesh().V()) 
                                    * runTime().deltaT().value() ;
                }
            } while(iterCtrl().loop()); //Check convergence of outer loop
            iterCtrl().write();
            ///////////////////- Nonlinear Iterations (Flux/Conductivity+Boundary Condition) converged //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (debug)
            {
                Info << "Timestep converged" << endl;
            }

            //Update fields
            pDot() = fvc::ddt(p());
            Info << "newPoroFluid evolved" << endl;

            return 0;
        }

        void newPoroFluid::writeFields(const Time &runTime)
        {
            Info << "writing PhysicalModel Fields" << endl;
            volScalarField p_(
                "p",
                p()());
            p_.write();

            volScalarField theta(
                "theta",
                poroHydraulic().S() * poroHydraulic().n());
            theta.write();

            volScalarField pHead(
                "pHead",
                p_ / poroHydraulic().magGamma());
            pHead.write();

            volScalarField kEff(
                "kEff",
                (poroHydraulic().kr() * fvc::reconstructMag(poroHydraulic().k()))());
            kEff.write();

            poroHydraulic().kr().write();

            volVectorField q(
                "q",
                fvc::reconstruct(phi()));
            q.write();

            poroFluidModel::writeFields(runTime);
        }

        poroHydraulicModel &newPoroFluid::poroHydraulic()
        const
        {
            if (!poroHydPtr_.valid())
            {
                FatalErrorIn("newPoroFluid::poroHydraulic()") << "poroHydraulicModel not initialized before first use" << endl;
            }

            return poroHydPtr_();
        }

        richardsLinearization &newPoroFluid::linearization()
        const
        {
            if (!richardsLinearizationPtr_.valid())
            {
               FatalErrorIn("newPoroFluid::linearization()") << "richardsLinearization not initialized before first use" << endl;
            }

            return richardsLinearizationPtr_();
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace poroFluidModels
} // End namespace Foam


// ************************************************************************* //
