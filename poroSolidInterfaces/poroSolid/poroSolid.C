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

#include "poroSolid.H"
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

        defineTypeNameAndDebug(poroSolid, 0);
        addToRunTimeSelectionTable(poroSolidInterface, poroSolid, dictionary);

        // * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //
        //- porosity change per unit time
        tmp<volScalarField> poroSolid::nDot(const volScalarField& b, const volVectorField& U)
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
        tmp<volScalarField> poroSolid::fixedStressStabil(const volScalarField& b, const volScalarField& impK)
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
        tmp<surfaceVectorField> poroSolid::q_relAcc(const surfaceScalarField& kf, const volVectorField& U)
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

        poroSolid::poroSolid(
            Time &runTime,
            const word &region)
            : poroSolidInterface(typeName, runTime, region),
              fixedStressStabil_(), // Pointer to stabilization term
              nDot_(), // Pointer to porosity change per unit time (exchange term)
              q_relAcc_(), // Pointer to rel. Acceleration term (exchange term)
              pRghSolidMesh_(), // p_rgh on solid Mesh (used in buoyancy less mode)
              pSolidMesh_() // p on solid Mesh (includes the buoyancy force on the solid)
        {
            // Sanity check if an appropriate poroSolid solver is used for the selected poroFluid solver
            if(poroFluid().name() == "varSatPoroFluid" || poroFluid().name() == "varSatPoroFluidHead")
            {
                Warning() << "'varSatPoroFluid/Head' should be used with 'varSatPoroSolid'" << endl;
            }
            
            // If we use the solid mesh for poroFluid calculations we only need to make the pressure available to the solid solver
            if(sharedMesh())
            {                
                // Checkin the pressure fields in the solid registry, so the porous material law wrapping class can find it later
                solidMesh().objectRegistry::checkIn(poroFluidRef().p());
                solidMesh().objectRegistry::checkIn(poroFluidRef().p_rgh());
            }
            else // For calculations where solid and fluid are on different meshes, we need to map the pressure field onto the solid mesh
            {
                // Map pressure fields onto solid field so the porous material law wrapping class can find it later
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

            }
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        bool poroSolid::evolve()
        {
            Info << "Preparing poroSolid solver" << endl;

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
                if(sharedMesh())
                {
                    // This is returned as "tmp" so we need to safe the data here otherwise it runs out of scope
                    const volScalarField impK(solid().mechanical().impK()); 
                    //calculate the exchange terms, the functions ndot, fixstressstabil and q_relAcc are defined above
                    nDot_.reset(nDot(b(),solid().U()).ptr());
                    fixedStressStabil_.reset(fixedStressStabil(b(),impK).ptr());
                    // TODO: Find solution to get this term working again
                    //q_relAcc_.reset(q_relAcc(poroFluid().poroHydraulic().kf()(),solid().U()).ptr());
                    // If we implicitly update the porosity, this is done now. 
                    if(!porosityConstantExplicit())
                    {
                        //Calculating n = n_start+div(D)
                        // Info: some of the compression might also be due to grain compression, 
                        // this is neglected here.
                        poroFluidRef().update_porosity(fvc::div(solid().D()),false); 
                    }
                }
                else // Different meshes for poroFluid and solid
                {
                    Info << "Mapping fields to poroFluid mesh";
                    // This is returned as "tmp" so we need to safe the data here otherwise it runs out of scope
                    const volScalarField impK(solid().mechanical().impK()); 
                    // Map the solid velocity onto fluid mesh
                    tmp<volVectorField> UFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().U());
                    // For the stabilization term, we also need the solid stiffness on the fluid mesh
                    tmp<volScalarField> tmpImpK(solidToPoroFluid().mapTgtToSrc(impK));

                    //calculate the exchange terms, the functions ndot, fixstressstabil and q_relAcc are defined above
                    nDot_.reset(nDot(b(),UFluidMesh()).ptr());
                    fixedStressStabil_.reset(fixedStressStabil(b(),tmpImpK()).ptr());
                    // TODO: Find solution to get this term working again
                    //q_relAcc_.reset(q_relAcc(poroFluid().poroHydraulic().kf()(),UFluidMesh()).ptr());

                    // The mapped velocity and stiffness are no longer needed, we can delete them to safe memory
                    tmpImpK.clear();
                    UFluidMesh.clear();

                    // If we implicitly update the porosity, this is done now. 
                    if(!porosityConstantExplicit())
                    {
                        //Maping the displacement field
                        tmp<volVectorField> DFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().D());
                        //Calculating n = n_start+div(D)
                        // Info: some of the compression might also be due to grain compression, 
                        // this is neglected here.
                        poroFluidRef().update_porosity(fvc::div(solid().D()),false); 
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
                    pRghSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluid().p_rgh())();
                    pSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluid().p());
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
                    // Info: some of the compression might also be due to grain compression, 
                    // this is neglected here.
                    poroFluidRef().update_porosity(fvc::div(solid().D()),false); 
                }
                else
                {
                    tmp<volVectorField> DFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().D());
                    //Calculating n = n_start+div(D)
                    // Info: some of the compression might also be due to grain compression, 
                    // this is neglected here.
                    poroFluidRef().update_porosity(fvc::div(solid().D()),false); 
                    DFluidMesh.clear();
                }
            }

            return true;
        }
        
            //- explicit coupling terms to pressure equation 
            const tmp<volScalarField> poroSolid::explicitCouplingDtoP() const
            {
                tmp<volScalarField> tSu(
                    new volScalarField(
                        nDot_() //- fvc::div(poroFluidMesh().Sf() & q_relAcc_())
                        )
                    );
                return tSu;
            }
            //- implicit coupling terms to pressure equation (get multiplied with pField)
            const tmp<volScalarField>  poroSolid::implicitCouplingDtoP() const
            {
                tmp<volScalarField> tSp(
                    new volScalarField(fixedStressStabil_())
                    );
                return tSp;
            }

        void poroSolid::writeFields(const Time &runTime)
        {
            poroFluidRef().writeFields(runTime);
            solidRef().writeFields(runTime);
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace poroSolidInteractions

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
