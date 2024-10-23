/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "abaqusUmatMohrCoulomb.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(abaqusUmatMohrCoulomb, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, abaqusUmatMohrCoulomb, linGeomMechLaw
    );

    // Declare fortran function prototypes
    extern "C"
    {
        // Note: all lowercase letters even if the fortran function has
        // uppercase letters
        void umat_(
            double STRESS[6],
            const double STATEV[], // double STATEV[6],
            double DDSDDE[6][6],
            const double *, // SSE,
            const double *, // SPD,
            const double *, // SCD,
            const double *, // RPL,
            const double *, // DDSDDT,
            const double *, // DRPLDE,
            const double *, // DRPLDT,
            const double STRAN[6],
            const double DSTRAN[6], // const double DSTRAN[6],
            const double *,         // TIME, (total time at beginning of increment)
            const double *,         // DTIME, (time increment)
            const double *,         // TEMP, (temperature at beginning of increment)
            const double *,         // DTEMP, (increment of temperature)
            const double *,         // PREDEF, (interp. values of predef. field variables)
            const double *,         // DPRED, (increment in predef. field variables)
            const double *,         // CMNAME, (no. of direct stress/strain components)
            const int *NDI,   // no. of direct stress/strain components
            const int *NSHR,  // no. of shear stress/strain components
            const int *NTENS,    //total no. of stress/strain components
            const int *NSTATV,      // no. of state variables
            const double PROPS[],   // array with material properties
            const int *NPROPS,    // (coordinates at point)
            const double *,      // COORDS,
            const double *,      // DROT, (coordinates at point)
            const double *,      // PNEWDT, (rotation increment vector)
            const double *,      // CELENT, (characteristic element length)
            const double *,      // DFGRD0, (deformation increment at start of step)
            const double *,      // DFGRD1, (deformation increment at end of step)
            const double *NOEL,  // NOEL, (element number)
            const double *NPT,   // NPT, (integration point number == 1)
            const double *,      // LAYER, (not used)
            const double *,      // KSPT, (not used)
            const double *KSTEP, // JSTEP,
            const double *KINC   // KINC
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::abaqusUmatMohrCoulomb::abaqusUmatMohrCoulomb(
    const word &name,
    const fvMesh &mesh,
    const dictionary &dict,
    const nonLinearGeometry::nonLinearType &nonLinGeom)
    : mechanicalLaw(name, mesh, dict, nonLinGeom),
      // rho_(dict.lookup("rho")),
      PROPS(),
      E_(dict.lookup("E")),
      nu_(dict.lookup("nu")),
      lambda_(
          planeStress()
              ? nu_ * E_ / ((1.0 + nu_) * (1.0 - nu_))
              : nu_ * E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_))),
      mu_(E_ / (2.0 * (1.0 + nu_))),
      K_(lambda_ + (2.0 / 3.0) * mu_),
      varPhi_(dict.lookup("frictionAngle")),
      c_(dict.lookup("cohesion")),
      varPsi_(dict.lookup("dilationAngle")),
      // stateVariables_(0),
      impK_(2 * mu_ + lambda_),
      DepsilonDevEq_(
          IOobject(
              "DEpsilonDevEq",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("0", dimless, 0)),
      activeYield_(
          IOobject(
              "activeYield",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("0", dimless, 0)),
      resSig_(0.0)
{
    // Force storage of strain old time
    epsilon().oldTime();

    PROPS[0] = E_.value();
    PROPS[1] = nu_.value();
    PROPS[2] = c_.value();
    PROPS[3] = varPhi_.value();
    PROPS[4] = varPsi_.value();

    // Initialise state varible fields

    // const scalarList stateVariablesInitialValues
    // (
    //     dict.lookup("stateVariablesInitialValues")
    // );

    // stateVariables_.setSize(stateVariablesInitialValues.size());

    // forAll(stateVariables_, fieldI)
    // {
    //     stateVariables_.set
    //     (
    //         fieldI,
    //         new volScalarField
    //         (
    //             IOobject
    //             (
    //                 "stateVariable" + Foam::name(fieldI + 1),
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::READ_IF_PRESENT,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh,
    //             dimensionedScalar
    //             (
    //                 "zero", dimless, stateVariablesInitialValues[fieldI]
    //             )
    //         )
    //     );

    //     // Force the old-time to be stored
    //     stateVariables_[fieldI].oldTime();
    // }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::abaqusUmatMohrCoulomb::~abaqusUmatMohrCoulomb()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::abaqusUmatMohrCoulomb::impK() const
{
    return tmp<volScalarField>(
        new volScalarField(
            IOobject(
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh(),
            impK_));
}

const Foam::dimensionedScalar &Foam::abaqusUmatMohrCoulomb::mu() const
{
    return mu_;
    ;
}

const Foam::dimensionedScalar &Foam::abaqusUmatMohrCoulomb::lambda() const
{
    return lambda_;
}

void Foam::abaqusUmatMohrCoulomb::correct(volSymmTensorField &sigma)
{

    sigma.storePrevIter();

    // Update epsilon
    updateEpsilon();

    // Update the increment of strain
    volSymmTensorField Depsilon("Depsilon",epsilon() - epsilon().oldTime());
    DepsilonDevEq_ = mag(dev(Depsilon));
    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn(
                "void Foam::abaqusUmatMohrCoulomb::"
                "correct(volSymmTensorField& sigma)")
                << "Not implemented for planeStress!" << abort(FatalError);
        }
    }

    // Initialise arrays to be passed to the fortran sub-routine, and then pass
    // them.
    // Note: many of the variables are not used by the abaqusUmatMohrCoulomb.C and
    // so we do not initialise them. In general, for other UMATS these may also
    // have to be initialised
    // BE CAREFUL: fortran expects column major indexing as opposed to row major
    // indexing, so we need to transpose all tensors before passing them and
    // receiving them
    // The following references are useful:
    // https://simplifiedfem.wordpress.com/about/tutorial-write-a-simple-umat-in-abaqus/
    // http://130.149.89.49:2080/v2016/books/sub/default.htm
    // http://130.149.89.49:2080/v2016/books/usb/default.htm?startat=pt01ch01s02aus02.html#usb-int-iconventions

    // Length of a stress tensor vector
    // Note: abaqus stores tensors as 1-D arrays (Voight notation_
    const int NTENS = 6;
    const int NDI = 3;    // number of direct components
    const int NSHR = 3;   // number of shear components
    const int NSTATV = 1; // stateVariablesI.size();

    // Material properties
    const int NPROPS = 5;

    // Internal field
    const symmTensorField &epsilonIOld = epsilon().oldTime().internalField();
    const symmTensorField &DepsilonI = Depsilon.internalField();
    const symmTensorField &sigmaIOld = sigma.oldTime().internalField();
    symmTensorField &sigmaI = sigma.internalFieldRef();

    forAll(epsilonIOld, cellI)
    {
        const symmTensor &sig = sigmaIOld[cellI];
        double STRESS[NTENS] =
            {
                sig.xx(), sig.yy(), sig.zz(),
                sig.xy(), sig.xz(), sig.yz()};
        double STATEV[NSTATV]; //Output only
        // forAll(stateVariables_, varI)
        // {
        //     STATEV[varI] = stateVariables_.internalField()[varI];
        // }
        double DDSDDE[NTENS][NTENS]; // the elastoplastic constitutive tensor (output)
        // double SSE[]; //  specific elastic strain energy
        // double SPD[]; // specific plastic dissipation
        // double SCD[]; // specific "creep" dissipation
        // double RPL[]; // volumetric heat generation 
        // double DDSDDT[NTENS]; //stress variation due to temperature*
        // double DRPLDE[NTENS]; // variation of RPL due to strain*
        // double DRPLDT[]; // variation of RPL due to temperature*
        const symmTensor &eps = epsilonIOld[cellI];
        double STRAN[NTENS] =
            {
                eps.xx(), eps.yy(), eps.zz(),
                2 * eps.xy(), 2 * eps.xz(), 2 * eps.yz()};
        const symmTensor &Deps = DepsilonI[cellI];
        double DSTRAN[NTENS] =
            {
                Deps.xx(), Deps.yy(), Deps.zz(),
                2 * Deps.xy(), 2 * Deps.xz(), 2 * Deps.yz()};
        // double TIME[2];
        // double DTIME;
        // double TEMP[];
        // double DTEMP[];
        // double PREDEF;
        // double DPRED;
        // double CMNAME[];
        // double COORDS[];
        // double DROT[];
        // double PNEWDT[];
        // double CELENT[];
        // double DFGRD0[];
        // double DFGRD1[];
        double NOEL = cellI;
        double NPT = 1.0;
        // double LAYER[];
        // double KSPT[];
        // double JSTEP[4] = 1;
        double KSTEP = 1;
        double KINC = 1;

        // Call Abaqus UMAT to calculate the stress
        double notImplemented = 0;
        umat_(
            STRESS,
            STATEV, // STATEV,
            DDSDDE,
            &notImplemented, // SSE,
            &notImplemented, // SPD,
            &notImplemented, // SCD,
            &notImplemented, // RPL,
            &notImplemented, // DDSDDT,
            &notImplemented, // DRPLDE,
            &notImplemented, // DRPLDT,
            STRAN,
            DSTRAN,          // DSTRAN,
            &notImplemented, // TIME,
            &notImplemented, // DTIME,
            &notImplemented, // TEMP,
            &notImplemented, // DTEMP,
            &notImplemented, // PREDEF,
            &notImplemented, // DPRED,
            &notImplemented, // CMNAME,
            &NDI,
            &NSHR,
            &NTENS,
            &NSTATV,
            PROPS,
            &NPROPS,
            &notImplemented, // COORDS,
            &notImplemented, // DROT,
            &notImplemented, // PNEWDT,
            &notImplemented, // CELENT,
            &notImplemented, // DFGRD0,
            &notImplemented, // DFGRD1,
            &NOEL,           // NOEL,
            &NPT,            // NPT,
            &notImplemented, // LAYER,
            &notImplemented, // KSPT,
            &KSTEP,          // JSTEP,
            &KINC            // KINC
        );

        // Retrieve the stress
        // If used, you would also retrieve the state variables
        sigmaI[cellI].xx() = STRESS[0];
        sigmaI[cellI].yy() = STRESS[1];
        sigmaI[cellI].zz() = STRESS[2];
        sigmaI[cellI].xy() = STRESS[3];
        sigmaI[cellI].xz() = STRESS[4];
        sigmaI[cellI].yz() = STRESS[5];
        activeYield_.internalFieldRef()[cellI] = STATEV[0];
    }

    // Boundary field
    forAll(epsilon().boundaryField(), patchI)
    {
        const symmTensorField &epsilonOldPatch = epsilon().oldTime().boundaryField()[patchI];
        const symmTensorField &DepsilonPatch = Depsilon.boundaryField()[patchI];
        const symmTensorField &sigmaOldPatch = sigma.oldTime().boundaryField()[patchI];
        symmTensorField &sigmaPatch = sigma.boundaryFieldRef()[patchI];

        forAll(epsilonOldPatch, faceI)
        {
            const symmTensor &sig = sigmaOldPatch[faceI];
            double STRESS[NTENS] =
                {
                    sig.xx(), sig.yy(), sig.zz(),
                    sig.xy(), sig.xz(), sig.yz()};
            double STATEV[NSTATV];
            // forAll(stateVariables_, varI)
            // {
            //     STATEV[varI] = stateVariables_.internalField()[varI];
            // }
            double DDSDDE[NTENS][NTENS];
            // double SSE[];
            // double SPD[];
            // double SCD[];
            // double RPL[];
            // double DDSDDT[NTENS];
            // double DRPLDE[NTENS];
            // double DRPLDT[];
            const symmTensor &eps = epsilonOldPatch[faceI];
            double STRAN[NTENS] =
                {
                    eps.xx(), eps.yy(), eps.zz(),
                    eps.xy(), eps.xz(), eps.yz()};
            const symmTensor &Deps = DepsilonPatch[faceI];
            double DSTRAN[NTENS] =
                {
                    Deps.xx(), Deps.yy(), Deps.zz(),
                    Deps.xy(), Deps.xz(), Deps.yz()};
            // double DSTRAN[NTENS];
            // double TIME[2];
            // double DTIME;
            // double TEMP[];
            // double DTEMP[];
            // double PREDEF;
            // double DPRED;
            // double CMNAME[];
            // double COORDS[];
            // double DROT[];
            // double PNEWDT[];
            // double CELENT[];
            // double DFGRD0[];
            // double DFGRD1[];
            double NOEL = faceI;
            double NPT = 1.0;
            // double LAYER[];
            // double KSPT[];
            // double JSTEP[4] = 1;
            double KSTEP = 1;
            double KINC = 1;

            // Call Abaqus UMAT to calculate the stress
            double notImplemented = 0;
            umat_(
                STRESS,
                STATEV, // STATEV,
                DDSDDE,
                &notImplemented, // SSE,
                &notImplemented, // SPD,
                &notImplemented, // SCD,
                &notImplemented, // RPL,
                &notImplemented, // DDSDDT,
                &notImplemented, // DRPLDE,
                &notImplemented, // DRPLDT,
                STRAN,
                DSTRAN,          // DSTRAN,
                &notImplemented, // TIME,
                &notImplemented, // DTIME,
                &notImplemented, // TEMP,
                &notImplemented, // DTEMP,
                &notImplemented, // PREDEF,
                &notImplemented, // DPRED,
                &notImplemented, // CMNAME,
                &NDI,
                &NSHR,
                &NTENS,
                &NSTATV,
                PROPS,
                &NPROPS,
                &notImplemented, // COORDS,
                &notImplemented, // DROT,
                &notImplemented, // PNEWDT,
                &notImplemented, // CELENT,
                &notImplemented, // DFGRD0,
                &notImplemented, // DFGRD1,
                &NOEL,           // NOEL,
                &NPT,            // NPT,
                &notImplemented, // LAYER,
                &notImplemented, // KSPT,
                &KSTEP,          // JSTEP,
                &KINC            // KINC
            );

            // Retrieve the stress
            // If used, you would also retrieve the state variables
            sigmaPatch[faceI].xx() = STRESS[0];
            sigmaPatch[faceI].yy() = STRESS[1];
            sigmaPatch[faceI].zz() = STRESS[2];
            sigmaPatch[faceI].xy() = STRESS[3];
            sigmaPatch[faceI].xz() = STRESS[4];
            sigmaPatch[faceI].yz() = STRESS[5];
            activeYield_.boundaryFieldRef()[patchI][faceI] = STATEV[0];
        }
    }
    // Calculate residual based on change in plastic strain increment
    resSig_ =
        gMax
        (
            mag
            (
                sigma.primitiveField()
                - sigma.prevIter().primitiveField()
            )/(SMALL + mag(sigma.primitiveField()-sigma.oldTime().primitiveField()))
        );
}

void Foam::abaqusUmatMohrCoulomb::correct(surfaceSymmTensorField &sigma)
{
    notImplemented("Foam::abaqusUmatMohrCoulomb::correct(surfaceSymmTensorField)");
}

Foam::scalar Foam::abaqusUmatMohrCoulomb::residual()
{
    return resSig_;
}

void Foam::abaqusUmatMohrCoulomb::updateTotalFields()
{
    
}


// ************************************************************************* //
