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

#include "varSatPoroSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicalModel.H"
#include "iterationControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace poroSolidInteractions
    {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(varSatPoroSolid, 0);
        addToRunTimeSelectionTable(poroSolidInterface, varSatPoroSolid, dictionary);

        // * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //
        //- porosity change per unit time
        tmp<volScalarField> varSatPoroSolid::nDot(const volScalarField& b, const volVectorField& U)
        {
            tmp<volScalarField> tnDot
            (
                new volScalarField
                (
                    "nDot",
                    b*fvc::div(U)
                )
            );
            return tnDot;
        }

        //- fixed stress stablization term
        tmp<volScalarField> varSatPoroSolid::fixedStressStabil(const volScalarField& b, const volScalarField& impK)
        {
            tmp<volScalarField> tStabil
            (
                new volScalarField
                (
                    "fixedStressSplitCoeff",
                    min
                    ( 
                        pow(b,2) 
                        / 
                        max(
                            impK,dimensionedScalar("0",impK.dimensions(),VGREAT)
                        ),
                        dimensionedScalar("0",dimless/impK.dimensions(),VSMALL)
                    )
                )
            );
            return tStabil;
        }

        //- fluxes arising from differencial acceleration (usually not significant)
        tmp<surfaceVectorField> varSatPoroSolid::q_relAcc(const surfaceScalarField& kf, const volVectorField& U)
        {
            tmp<surfaceVectorField> tq
            (
                new surfaceVectorField
                (
                    "q_relAcc",
                    kf * fvc::interpolate(fvc::ddt(U)/mag(poroFluid().g()))
                )
            );
            return tq;
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
        varSatPoroSolid::varSatPoroSolid(
            Time &runTime,
            const word &region)
            : poroSolidInterface(typeName, runTime, region),
              fixedStressStabil_(), // Pointer to stabilization term
              nDot_(), // Pointer to porosity change per unit time (exchange term)
              q_relAcc_(), // Pointer to rel. Acceleration term (exchange term)
              pRghSolidMesh_(), // p_rgh on solid Mesh (used in buoyancy less mode)
              pSolidMesh_(), // p on solid Mesh (includes the buoyancy force on the solid)
              SSolidMesh_(), // saturation on solid Mesh 
              nSolidMesh_(), // porosity on solid Mesh 
              magGammaW_(dimensionedScalar("gammaw",dimless,1.0)) // this to facilitate pressureHead and p_rgh based solvers in one class
        {
            // Sanity check if an appropriate poroSolid solver is used for the selected poroFluid solver
            if(poroFluid().name() == "poroFluid")
            {
                FatalErrorInFunction() << "'poroFluid' is not suitable for variably saturated calculations" << endl;
            }

            // this to facilitate pressureHead and p_rgh based solvers in one class
            // we will use either 1 or specific weight of water to bring all fields to the chosen pressure units
            if(pField().dimensions()==dimLength)
            {
                const dimensionedVector& gammaW = poroFluid().lookupObject<UniformDimensionedField<vector>>("gamma_water");
                magGammaW_ = mag(gammaW).value();
                magGammaW_.dimensions() = gammaW.dimensions();
            }
            else if(pField().dimensions()!=dimPressure)
            {
                FatalError() << "pore pressure field is neither in head nor in pressure dimensions!"
                             << endl;
            }

            // Get reference to saturation
            const volScalarField& S = poroFluid().lookupObject<volScalarField>("S");

            // If we use the solid mesh for poroFluid calculations we only need to make the pressure available to the solid solver
            if(sharedMesh())
            {
                // Checkin the pressure fields in the solid registry, so the porous material law wrapping class can find it later
                solidMesh().objectRegistry::checkIn(poroFluidRef().p().ref());
                solidMesh().objectRegistry::checkIn(poroFluidRef().p_rgh().ref());
                solidMesh().objectRegistry::checkIn(const_cast<volScalarField&>(S));
                solidMesh().objectRegistry::checkIn(const_cast<volScalarField&>(poroFluidRef().n()));
            }
            else
            {

                // Map pressure fields onto solid field so the porous material law wrapping class can find it later
                // This is done now, and the fields are kept in memory so that other classes can use them at any time
                // e.g. the material model needs them very early bc they are used to calculate total density of solid+water mixture
                pSolidMesh_.reset(
                    new volScalarField(
                        IOobject(
                            "p",
                            runTime.timeName(),
                            solidMesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE),
                            solidToPoroFluid().mapSrcToTgt(poroFluid().p())()));

                pRghSolidMesh_.reset(
                    new volScalarField(
                        IOobject(
                            "p_rgh",
                            runTime.timeName(),
                            solidMesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE),
                            solidToPoroFluid().mapSrcToTgt(poroFluid().p_rgh())()));

                SSolidMesh_.reset(
                    new volScalarField(
                        IOobject(
                            "S",
                            runTime.timeName(),
                            solidMesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE),
                            solidToPoroFluid().mapSrcToTgt(S)(),
                            "zeroGradient"));

                nSolidMesh_.reset(
                    new volScalarField(
                        IOobject(
                            "n",
                            runTime.timeName(),
                            solidMesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE),
                            solidToPoroFluid().mapSrcToTgt(poroFluidRef().n())(),
                            "zeroGradient"));
            }
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        bool varSatPoroSolid::evolve()
        {
            Info << "Preparing varSatPoroSolid solver" << endl;

            SolverPerformance<scalar> solverPerfp;
            SolverPerformance<vector>::debug = 0;

            // Reset the residuals from last timestep
            couplingControl().reset();

            //#######- Pressure-displacement coupling outer loop #########
            do
            {
                //## First initialize of coupling terms for the fluid

                // For shared mesh calculations (solid and poroFluid use the same mesh) 
                // We dont need any mesh to mesh mapping, we can simply initilize the fields on the solid mesh
                // Get reference to saturation
                const volScalarField& SFluidMesh = poroFluid().lookupObject<volScalarField>("S");
                if(sharedMesh())
                {
                    // This is returned as "tmp" so we need to safe the data here otherwise it runs out of scope
                    const volScalarField impK(solid().mechanical().impK()/magGammaW_);
                    // We combine bishop parameter (effect of partial saturation on the coupling)
                    // and biot parameter (effect of compressible solid consituents)
                    tmp<volScalarField> bishopBiot(SFluidMesh*b());
                    //calculate the exchange terms, the functions ndot, fixstressstabil and q_relAcc are defined above
                    nDot_.reset(nDot(bishopBiot(),solid().U()).ptr());
                    fixedStressStabil_.reset(fixedStressStabil(b(),impK).ptr());
                    // TODO: Find solution to get this term working again
                    //q_relAcc_.reset(q_relAcc(poroHydraulic().kEfff()(),solid().U()).ptr());
                    // If we implicitly update the porosity, this is done now. 
                    if(!porosityConstantExplicit())
                    {
                        //Calculating n = n_start+div(D)
                        // Info: some of the compression might also be due to grain compression, 
                        // this is neglected here.
                        poroFluidRef().update_porosity(fvc::div(solid().D()),false); 
                    }
                }
                else
                {
                    Info << "Mapping fields to poroFluid mesh" << endl;
                    // This is returned as "tmp" so we need to safe the data here otherwise it runs out of scope
                    const volScalarField impK(solid().mechanical().impK()/magGammaW_);
                    // Map the solid velocity onto fluid mesh
                    tmp<volVectorField> UFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().U());
                    // For the stabilization term, we also need the solid stiffness on the fluid mesh
                    tmp<volScalarField> tmpImpK(solidToPoroFluid().mapTgtToSrc(impK)());
                    // We combine bishop parameter (effect of partial saturation on the coupling)
                    // and biot parameter (effect of compressible solid consituents)
                    tmp<volScalarField> bishopBiot(SFluidMesh*b());

                    //calculate the exchange terms, the functions ndot, fixstressstabil and q_relAcc are defined above
                    nDot_.reset(nDot(bishopBiot(),UFluidMesh()).ptr());
                    fixedStressStabil_.reset(fixedStressStabil(b(),tmpImpK).ptr());
                    // TODO: Find solution to get this term working again
                    //q_relAcc_.reset(q_relAcc(poroHydraulic().kEfff()(),UFluidMesh()).ptr());
                    
                    // The mapped velocity and stiffness are no longer needed, we can delete them to safe memory
                    tmpImpK.clear();
                    UFluidMesh.clear();
                    bishopBiot.clear();

                    // If we implicitly update the porosity, this is done now. 
                    if(!porosityConstantExplicit())
                    {
                        //Maping the displacement field
                        tmp<volVectorField> DFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().D());
                        //Calculating n = n_start+div(D) 
                        // Info: some of the compression might also be due to grain compression, 
                        // this is neglected here.
                        poroFluidRef().update_porosity(fvc::div(DFluidMesh),false); 
                        nSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluidRef().n())();
                        // Clear mapped D to safe memory
                        DFluidMesh.clear();
                    }
                }

                //- Evolving the fluid solver
                poroFluidRef().evolve();

                // Delete mechanic to hydraulic coupling terms to safe memory
                nDot_.clear();
                fixedStressStabil_.clear();
                q_relAcc_.clear();

                //- Preparing the hydraulic to mechanic source terms for coupling
                // If we share the mesh for solid and poroFluid, the pressures are already checked in
                // The poroMechanicalLaw class will take care of the coupling terms
                // on different meshes, we first need to map the pressure fields from poroFluid to solid mesh
                if (!sharedMesh())
                {
                    // Get reference to saturation
                    const volScalarField& S = poroFluid().lookupObject<volScalarField>("S");
                    pRghSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluid().p_rgh())();
                    pSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluid().p());
                    SSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(S);
                }

                // To update the solid+water total density that is stored in the solid material law
                // we need to envoke this funciton
                // This is only done from the second iteration onwards, we want the inital to be the provided
                if(couplingControl().index()>0)
                {
                    solidRef().recalculateRho();
                }

                //- Evolving solid solver
                solidRef().evolve();

            } while (couplingControl().loop());

            // If we wont to write out iteration metrics, this function will do it
            couplingControl().write();

            Info << "Coupling Evolved" << endl;

            // If porosity is changing we do this here.
            // This is the only time we update porosity in the timestep for explicit porosity
            // and we do it here again for implicit porosity to make sure its uptodate
            if(!porosityConstant())
            {
                if(sharedMesh())
                {
                    //Calculating n = n_start+div(D)
                    poroFluidRef().update_porosity(fvc::div(solid().D()),false); 
                }
                else
                {
                    tmp<volVectorField> DFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().D());
                    //Calculating n = n_start+div(D)
                    poroFluidRef().update_porosity(fvc::div(DFluidMesh),false); 
                    nSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluidRef().n())();
                    DFluidMesh.clear();
                }
            }

            return true;
        }

        //- explicit coupling terms to pressure equation
        const tmp<volScalarField> varSatPoroSolid::explicitCouplingDtoP() const
        {
            tmp<volScalarField> tSu(
                new volScalarField(
                    nDot_() //- fvc::div(poroFluidMesh().Sf() & q_relAcc_())
                    )
                );
            return tSu;
        }
        //- implicit coupling terms to pressure equation (get multiplied with pField)
        const tmp<volScalarField>  varSatPoroSolid::implicitCouplingDtoP() const
        {
            tmp<volScalarField> tSp(
                new volScalarField(fixedStressStabil_())
                );
            return tSp;
        }

        void varSatPoroSolid::writeFields(const Time &runTime)
        {
            solid().rho().write();

            poroFluidRef().writeFields(runTime);
            solidRef().writeFields(runTime);
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace poroSolidInteractions

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
