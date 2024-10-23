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

#include "varSatPoroFluid.H"
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

        defineTypeNameAndDebug(varSatPoroFluid, 0);
        addToRunTimeSelectionTable(poroFluidModel, varSatPoroFluid, dictionary);

        // * * * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

        void varSatPoroFluid::makeKEffbyGammaf()
        {
            kEffbyGammafPtr_.reset(
                new surfaceScalarField(
                IOobject(
                      "kEffbyGammaf",
                      runTime().timeName(),
                      p_rgh().db(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  p_rgh().mesh(),
                  dimensionedScalar("", dimensionSet(-1, 3, 1, 0, 0, 0, 0), 0.0)
            ));
        }
    
        void varSatPoroFluid::makePoroHydraulic()
        {
            poroHydPtr_.reset(new varSatPoroHydraulicModel(p_rgh(), g()));
        }

        // * * * * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * * * //

        surfaceScalarField& varSatPoroFluid::kEffbyGammaf()
        {
            if (!kEffbyGammafPtr_.valid())
            {
                makeKEffbyGammaf();
            }

            return kEffbyGammafPtr_.ref();
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        varSatPoroFluid::varSatPoroFluid(
            Time &runTime,
            const word &region,
            const bool sharedmesh)
            : variablySaturatedPoroFluid(typeName, runTime, "p_rgh", region, sharedmesh),
              Ss_ //- Storage coefficient
              (
                (
                    IOobject
                    (
                        "Ss", 
                        runTime.timeName(),
                        p_rgh().db(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    poroHydraulic().Ss(n(),p_rgh())
                )
              ),
              kEffbyGammafPtr_(), // ... effective hydraulic conductivity (saturated * unsaturated)/gamma_w
              steadyState_(word(mesh().ddtScheme("ddt(C,p_rgh)")) == "steadyState")
        {
            pRGHisRequired();
            //- Check if we have a steady state case (in which case underrelaxation for kr can be activated, otherwise not)
            if (!steadyState_ && mesh().relaxField("kEffbyGammaf") && mesh().fieldRelaxationFactor("kEffbyGammaf") != 1.0)
            {
                WarningInFunction() << " k should only be relaxed for steady-state calculations!!!"
                //FatalErrorIn("varSatPoroFluidHead::varSatPoroFluidHead") << " k should only be relaxed for steady-state calculations!!!"
                                                                   << endl;
            }

            // Update total pore pressure
            p() = p_rgh() + poroHydraulic().p_Hyd();

            // Now that p is initalized we can initalize these object:
            SPtr_.reset(poroHydraulic().S(p()).ptr());

            // Initialize iteration control
            iterCtrl();
            // Check if iteration control checks mass balance,
            // if yes, we need to provide this field
            checkMassBalance();

            if (debug)
            {
                Info << p_rgh().db().sortedNames()
                     << endl;
                    
                n().write();
                poroHydraulic().k()().write();
                Ss_.write();
                S().write();
                linearization().C().write();
            }

            Info << "Done generating variable saturated poroFluidModel (p_rgh)" << endl;
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        bool varSatPoroFluid::evolve()
        {

            // Initialize performance tracking
            // and silence automatic solver output

            SolverPerformance<scalar> solverPerfp;
            SolverPerformance<scalar>::debug = 0;

            // Resetting all residuals
            iterCtrl().reset();

            Info << "Evolving fluid model: " << this->type() << endl;

            ///////////- Save Fields from previous timeStep ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Update total pore pressure
            p() = p_rgh() + poroHydraulic().p_Hyd();
            // Update saturation
            S() = poroHydraulic().S(p());
            // Update the effective hydraulic conductivity
            kEffbyGammaf() = poroHydraulic().kEfff(p()) / poroHydraulic().magGamma();
            // Check if poroHydraulicModel changes storage
            // If yes, update storage
            if(poroHydraulic().updatesSs())
            {
                Ss_ = poroHydraulic().Ss(n(),p());
            }

            ///////////////////- Update Flux ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Update hydraulic gradient for use in boundary conditions and models
            hydraulicGradient() = fvc::grad(p_rgh()) / poroHydraulic().magGamma();
            // Update volumetric flux for use in boundary conditions
            phi() = -kEffbyGammaf() * (mesh().magSf() * fvc::snGrad(p_rgh()));


            ///////////- Nonlinear Iterations (Flux/Conductivity+Boundary Condition) ////////////////////////////////////////////////////////////////////////////////////////////////////
            do
            {

                // Clear previous outer iterations convergence data
                const_cast<dictionary &>(mesh().solverPerformanceDict()).clear();

                // Store prev fields for relaxation
                p_rgh().storePrevIter();

                ///////////////////- Initialize/Update the coefficients //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                linearization().initalize(p(), p_rgh());
                // In case Casulli's Scheme is used, p() is capped initially,
                // this needs to be transfered to p_rgh
                // p_rgh().correctBoundaryConditions();
                p_rgh() = p() - poroHydraulic().p_Hyd(); 

                //////////////////- Nonlinear Iterations (Saturation)  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                do
                {

                    // Debug statements to see if coupling works
                    if (debug)
                    {
                        Info << "Assembling pEqn"
                             << endl;

                        for (fv::option &source : static_cast<fv::optionList &>(fvOptions()))
                        {
                            Info << "Source Terms " << source.name() << " is active? " << source.active() << endl;
                        }
                        SolverPerformance<scalar>::debug = debug;
                    }

                    // Assemble pressure equation
                    fvScalarMatrix pEqn(
                        n() * linearization().ddtS(S(), p_rgh()) //- Change in Saturation
                            + S() * fvm::ddt(Ss_, p_rgh())      //- Storage (e.g. compressibility of fluid)
                            -fvm::laplacian(kEffbyGammaf(), p_rgh())           //- Implicit darcian fluxes
                        ==
                            fvOptions()(Ss_, p_rgh()) //- Source terms (e.g. volumetric deformation, wells etc.)
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

                    // Solve System
                    solverPerfp = pEqn.solve();

                    fvOptions().correct(p_rgh());
                    p() = p_rgh() + poroHydraulic().p_Hyd();// Update the total pressure with new excess pressure

                    phi() = pEqn.flux();

                    if (debug)
                    {
                        Info << "Solved pEqn" << endl;
                    }

                } while (!linearization().checkConvergedAndUpdate(p(), p_rgh())); // Check Convergence and Update Internal Veriables
                //////////////////- Nonlinear Iterations (Saturation) converged  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Relax p_rgh
                p_rgh().relax();
                p() = p_rgh() + poroHydraulic().p_Hyd();
                hydraulicGradient() = fvc::grad(p_rgh()) / poroHydraulic().magGamma();

                //////////////////- Output first and latest linear solver performance  ////////////////////////////////////////////////////////////////////////////////////////////////
                autoPtr<List<SolverPerformance<scalar>>> sp(new List<SolverPerformance<scalar>>(mesh().solverPerformanceDict().findEntry("p_rgh")->stream()));
                if (sp().size() == 1)
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
                
                // Update Saturation
                S() = poroHydraulic().S(p());

                // If this is a steady state calculation, dont relax kEff (leads to wrong solutions)
                // In any way update kEff
                if (steadyState_)
                {
                    kEffbyGammaf() .storePrevIter();
                    kEffbyGammaf() = poroHydraulic().kEfff(p()) / poroHydraulic().magGamma();
                    kEffbyGammaf().relax();
                }
                else
                {
                    kEffbyGammaf() = poroHydraulic().kEfff(p()) / poroHydraulic().magGamma();
                }

                // If saturationLaw updates Ss, then do it now
                if(poroHydraulic().updatesSs())
                {
                    Ss_ = poroHydraulic().Ss(n(),p());
                }

                p_rgh().correctBoundaryConditions();

            } while (iterCtrl().loop()); // Check convergence of outer loop

            iterCtrl().write();
            ///////////////////- Nonlinear Iterations (Flux/Conductivity+Boundary Condition) converged //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (debug)
            {
                Info << "Timestep converged" << endl;
            }

            // Update fields
            pDot() = fvc::ddt(p());
            Info << "varSatPoroFluid evolved" << endl;

            return 0;
        }

        scalar varSatPoroFluid::newDeltaT()
        {
            scalar maxCv = max(poroHydraulic().kEff(p()) / ((Ss_ + linearization().C())*poroHydraulic().magGamma())).value();
            Info << "max c_v: " << maxCv << endl;
            scalar minDeltaX = min(mesh().V() / fvc::surfaceSum(mag(mesh().Sf())));
            scalar wishedTimeStep = CoNumber_ * pow(minDeltaX, 2) / maxCv;
            scalar maxDeltaTFact = wishedTimeStep / runTime().deltaTValue();
            scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);
            return deltaTFact * runTime().deltaTValue();
        }

        void varSatPoroFluid::writeFields(const Time &runTime)
        {
            Info << "writing PhysicalModel Fields" << endl;

            volScalarField theta
            (
                "theta",
                S() * n()
            );
            theta.write();

            volScalarField pHead
            (
                "pHead",
                p() / poroHydraulic().magGamma()
            );
            pHead.write();

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

            volVectorField q
            (
                "q",
                fvc::reconstruct(phi())
            );
            q.write();

            poroFluidModel::writeFields(runTime);
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace varSatPoroFluidModels
} // End namespace Foam

// ************************************************************************* //
