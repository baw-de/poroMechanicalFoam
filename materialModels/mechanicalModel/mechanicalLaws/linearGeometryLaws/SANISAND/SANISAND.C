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

    You should have received a_ copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SANISAND.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SANISAND, 0);
    addToRunTimeSelectionTable(
        mechanicalLaw, SANISAND, linGeomMechLaw);
}

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

void Foam::SANISAND::calculateStress(
    symmTensor &stress,
    const scalar &dEpsV,
    const symmTensor &dEpsD,
    const scalar &poreRatio,
    const scalar &POld,
    const symmTensor &sOld,
    const symmTensor &alphaOld,
    const symmTensor &FOld,
    const scalar &mOld,
    scalar &dEpsPV,
    symmTensor &dEpsPD,
    scalar &PNew,
    symmTensor &sNew,
    symmTensor &alphaParam,
    scalar &mParam,
    symmTensor &FParam,
    scalar &KParam,
    scalar &GParam,
    scalar &psiParam,
    scalar &yieldFlag) const
{
    dEpsPV = 0.0;
    dEpsPD = symmTensor::zero;
    yieldFlag = 0;
    alphaParam = alphaOld;
    mParam = mOld;
    FParam = FOld;

    int i(0);
    scalar f = GREAT;

    do
    {
        if (debug == 2)
        {
            Info << "Calculate new P_, s_, psi_" << endl;
            Info << "Increment of epsilon volumetric: " << dEpsV << nl
                 << "Increment of epsilon plastic volumetric: " << dEpsPV << endl;
            Info << "Old P_: " << POld << endl;
            Info << "Old s_: " << sOld << endl;
        }

        if (debug == 3)
        {
            if(POld < SMALL)
            {
                WarningIn("Foam::SANISAND::calculateStress()") << "Mean effective pressure is zero or negative: " << POld << endl
                                                << "This is not allowed in Sanisand model!!" << endl;
            }
        }
        scalar PIntermediate = Foam::pow(
            POld, 1 - a_.value());
        scalar IntermediateTerms = (K0_.value() / Foam::pow(Pref_.value(), a_.value()) * (1.0 - a_.value()) * -(dEpsV - dEpsPV));
        PNew = Foam::pow(
            PIntermediate + IntermediateTerms,
            1.0 / (1.0 - a_.value()));
        sNew = sOld + 2.0 * GParam * -(dEpsD - dEpsPD);
        psiParam = poreRatio - (e_cref_.value() - LAMBDA_.value() * Foam::log(PNew / Pref_.value()));

                
        if (debug == 2)
        {
            Info << "New P_: " << PNew << endl;
            Info << "New s_: " << sNew << endl;
            Info << "New psi_: " << psiParam << endl;
        }

        f = mag(sNew - PNew * alphaParam) - Foam::sqrt(2.0 / 3.0) * mParam * PNew;

        if (f > 0 && (mag(f / Pref_.value()) > tol_f_))
        {

            symmTensor r_ = sNew / PNew - alphaParam;

            scalar J_ = Foam::sqrt(0.5 * tr(r_ & r_));
            scalar S_(0.0);
            if (tr(r_ & r_ & r_) < 0)
            {
                S_ = -Foam::pow(-1.0 / 3.0 * tr(r_ & r_ & r_), 1.0 / 3.0);
            }
            else
            {
                S_ = Foam::pow(1.0 / 3.0 * tr(r_ & r_ & r_), 1.0 / 3.0);
            }

            scalar cos3theta = 3.0 * Foam::sqrt(3.0) / 2.0 * (S_ / J_) * (S_ / J_) * (S_ / J_);
            scalar gc = (2.0 * c_ / ((1.0 + c_) - (1.0 - c_) * cos3theta)).value();
            scalar gb = (2.0 * c_b_ / ((1.0 + c_b_) - (1.0 - c_b_) * cos3theta)).value();
            scalar gd = (2.0 * c_d_ / ((1.0 + c_d_) - (1.0 - c_d_) * cos3theta)).value();

            symmTensor n_ = (sNew - PNew * alphaParam) / mag(sNew - PNew * alphaParam);

            symmTensor alpha_c = Foam::sqrt(2.0 / 3.0) * (gc * Mc_.value() - mParam) * n_;
            symmTensor alpha_b = Foam::sqrt(2.0 / 3.0) * (gc * Mc_.value() + gb * k_cb_.value() * (psiParam <= 0 ? (-psiParam) : 0) - mParam) * n_;
            symmTensor alpha_d = Foam::sqrt(2.0 / 3.0) * (gc * Mc_.value() + gd * k_cd_.value() * psiParam - mParam) * n_;

            symmTensor b = alpha_b - alphaParam;
            symmTensor d = alpha_d - alphaParam;

            scalar N = (alphaParam && n_) + 2.0 / 3.0 * mParam;
            scalar A = A0_.value() * (1.0 + ((FParam && n_) >= 0 ? (FParam && n_) : 0));

            scalar D = A * (d && n_);

            KParam = K0_.value() * Foam::pow(PNew / Pref_.value(), a_.value());
            GParam = G0_.value() * Foam::pow(PNew / Pref_.value(), a_.value());

            scalar bref = 2.0 * Foam::sqrt(2.0 / 3.0) * (gc * Mc_.value() + k_cb_.value() * (psiParam <= 0 ? (-psiParam) : 0) - mParam);
            scalar h = h0_.value() * (mag(b && n_)) / (bref - mag(b && n_));

            scalar Kp = PNew * (h * (b && n_) + Foam::sqrt(2.0 / 3.0) * cm_.value() * (1.0 + e0_.value()) * D);
            scalar dlambda = f / (2.0 * GParam - N * D * KParam + Kp);

            dEpsPV -= plast_relax_ * dlambda * D;
            dEpsPD -= plast_relax_ * dlambda * n_;

            alphaParam += dlambda * h * b;
            mParam += dlambda * cm_.value() * (1.0 + e0_.value()) * D;
            FParam -= dlambda * Cf_.value() * (D <= 0 ? (-D) : 0) * (Fmax_.value() * n_ + FParam);
            yieldFlag = 1;

            i++;
        }
    } while ((f > 0 && mag(f / Pref_.value()) > tol_f_) && (i < nCorr_));

    stress = -(PNew * symmTensor::I + sNew);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::SANISAND::SANISAND(
    const word &name,
    const fvMesh &mesh,
    const dictionary &dict,
    const nonLinearGeometry::nonLinearType &nonLinGeom)
    : mechanicalLaw(name, mesh, dict, nonLinGeom),
      rho_(dict.lookup("rho")),
      e0_(dict.lookup("e0")),
      // Hypoelasticity
      K0_(dict.lookup("K0")),
      G0_(dict.lookup("G0")),
      lambda0_(K0_ - 2 / 3 * G0_),
      e_cref_(dict.lookup("e_cref")),
      Pref_(dict.lookup("Pref")),
      a_(dict.lookup("a")),
      // Critical state line
      Mc_(dict.lookup("Mc")),
      Me_(dict.lookup("Me")),
      LAMBDA_(dict.lookup("lambda")),
      // hardening
      A0_(dict.lookup("A0")),
      Cf_(dict.lookup("Cf")),
      Fmax_(dict.lookup("Fmax")),
      h0_(dict.lookup("h0")),
      m0_(dict.lookup("m0")),
      cm_(dict.lookup("cm")),
      k_cb_(dict.lookup("k_cb")),
      k_eb_(dict.lookup("k_eb")),
      k_cd_(dict.lookup("k_cd")),
      k_ed_(dict.lookup("k_ed")),
      c_("c", dimless, (Me_ / Mc_).value()),
      c_b_("c_b", dimless, (k_eb_ / k_cb_).value()),
      c_d_("c_d", dimless, (k_ed_ / k_cd_).value()),
      tol_f_(readScalar(dict.lookup("toleranceYield"))),
      plast_relax_(readScalar(dict.lookup("relaxationPlasticity"))),
      nCorr_(readScalar(dict.lookup("nCorr"))),
#include "createFields.H"
      activeYield_(
          IOobject(
              "activeYield",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("0", dimless, 0))
{
    //- Force Field to safe oldTime
    e_.storeOldTime();
    const volSymmTensorField &sigma = mesh.lookupObject<volSymmTensorField>("sigma");
    if (debug == 1)
    {
        Info << "Generating SANISAND stress fields" << endl;
    }
    // P_.oldTime() = -(1.0 / 3.0) * tr(sigma.oldTime());
    P_ = -(1.0 / 3.0) * tr(sigma);
    if(min(P_).value() < SMALL)
    {
        WarningIn("Foam::SANISAND::SANISAND()") << "Mean effective pressure is zero or negative: " << min(P_) << endl
                                                << "This is not allowed in Sanisand model!!" << endl; 
    }
    // s_.oldTime() = -Foam::dev(sigma.oldTime());
    s_ = -Foam::dev(sigma);
    if (debug == 1)
    {
        Info << "Generating SANISAND stiffness fields" << endl;
        if (debug == 2)
        {
            Info << P_ << endl;
        }
    }
    K_ = K0_ * Foam::pow(P_ / Pref_, a_);
    mu_ = G0_ * Foam::pow(P_ / Pref_, a_);
}

/*  lambda_
  (
      planeStress()
    ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
    : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
  ),*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SANISAND::~SANISAND()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::SANISAND::rho() const
{
    tmp<volScalarField> tresult(
        new volScalarField(
            IOobject(
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName));

    tresult.ref().correctBoundaryConditions();

    return tresult;
}

Foam::tmp<Foam::volScalarField>
Foam::SANISAND::impK() const
{
    return tmp<volScalarField>(
        new volScalarField(
            IOobject(
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            4 / 3 * mu_ + K_));
}

const Foam::dimensionedScalar &Foam::SANISAND::mu() const
{
    return G0_;
}

const Foam::dimensionedScalar &
Foam::SANISAND::lambda() const
{
    return lambda0_;
}

void Foam::SANISAND::correct(volSymmTensorField &sigma)
{
    if (debug == 1)
    {
        Info << "Start calculating new sigma" << endl;
    }
    // Store the previous iteration of sigma to allow under-relaxation and also
    // to calculate the material residual
    sigma.storePrevIter();
    DSigma_.storePrevIter();

    // Calculate the increment of total strain
    if (debug == 1)
    {
        Info << "calculating DEpsilon_" << endl;
    }

    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField &gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        DEpsilon_ = symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField &gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate gradient of displacement increment
        const volTensorField gradDD = gradD - gradD.oldTime();

        DEpsilon_ = symm(gradDD);
    }

    DEpsilonV_ = tr(DEpsilon_);
    DEpsilonD_ = dev(DEpsilon_);

    if (debug == 1)
    {
        Info << "Calculating new pore-Ratio" << endl;
    }
    e_ = e_.oldTime() + (1 + e0_) * (DEpsilonV_);

    // Take a_ reference to internal fields for efficiency
    if (debug == 1)
    {
        Info << "Take a_ reference to internal fields for efficiency" << endl;
    }
    symmTensorField &sigmaI = sigma.primitiveFieldRef();
    scalarField &activeYieldI = activeYield_.primitiveFieldRef();
    const scalarField &dEpsV_I = DEpsilonV_.primitiveField();
    const symmTensorField &dEpsD_I = DEpsilonD_.primitiveField();
    const scalarField &eI = e_.primitiveField();
    const scalarField &POldI = P_.oldTime().primitiveField();
    const symmTensorField &sOldI = s_.oldTime().primitiveField();
    const symmTensorField &alphaOldI = alpha_.oldTime().primitiveField();
    const symmTensorField &FOldI = F_.oldTime().primitiveField();
    const scalarField &mOldI = m_.oldTime().primitiveField();
    scalarField &dEpsPVI = DEpsilonPV_.primitiveFieldRef();
    symmTensorField &dEpsPDI = DEpsilonPD_.primitiveFieldRef();
    scalarField &PI = P_.primitiveFieldRef();
    symmTensorField &sI = s_.primitiveFieldRef();
    symmTensorField &alphaI = alpha_.primitiveFieldRef();
    scalarField &mI = m_.primitiveFieldRef();
    symmTensorField &FI = F_.primitiveFieldRef();
    scalarField &KI = K_.primitiveFieldRef();
    scalarField &GI = mu_.primitiveFieldRef();
    scalarField &psiI = psi_.primitiveFieldRef();

    // Correct sigma internal field

    if (debug == 1)
    {
        Info << "calculating stress for internal field" << endl;
    }
    forAll(sigmaI, cellI)
    {
        calculateStress(sigmaI[cellI],
                        dEpsV_I[cellI],
                        dEpsD_I[cellI],
                        eI[cellI],
                        POldI[cellI],
                        sOldI[cellI],
                        alphaOldI[cellI],
                        FOldI[cellI],
                        mOldI[cellI],
                        dEpsPVI[cellI],
                        dEpsPDI[cellI],
                        PI[cellI],
                        sI[cellI],
                        alphaI[cellI],
                        mI[cellI],
                        FI[cellI],
                        KI[cellI],
                        GI[cellI],
                        psiI[cellI],
                        activeYieldI[cellI]);
    }

    // Correct sigma on the boundary patches
    if (debug == 1)
    {
        Info << "Correct sigma on the boundary patches" << endl;
    }
    forAll(sigma.boundaryField(), patchI)
    {
        // Take references to the boundary patches for efficiency
        if (debug == 1)
        {
            Info << "Take references to the boundary patches for efficiency" << endl;
        }

        symmTensorField &sigmaP = sigma.boundaryFieldRef()[patchI];
        scalarField &activeYieldP = activeYield_.boundaryFieldRef()[patchI];
        const scalarField &dEpsV_Patch = DEpsilonV_.boundaryField()[patchI];
        const symmTensorField &dEpsD_Patch = DEpsilonD_.boundaryField()[patchI];
        const scalarField &ePatch = e_.boundaryField()[patchI];
        const scalarField &POldPatch = P_.oldTime().boundaryField()[patchI];
        const symmTensorField &sOldPatch = s_.oldTime().boundaryField()[patchI];
        const symmTensorField &alphaOldPatch = alpha_.oldTime().boundaryField()[patchI];
        const symmTensorField &FOldPatch = F_.oldTime().boundaryField()[patchI];
        const scalarField &mOldPatch = m_.oldTime().boundaryField()[patchI];
        scalarField &dEpsPVPatch = DEpsilonPV_.boundaryFieldRef()[patchI];
        symmTensorField &dEpsPDPatch = DEpsilonPD_.boundaryFieldRef()[patchI];
        scalarField &PPatch = P_.boundaryFieldRef()[patchI];
        symmTensorField &sPatch = s_.boundaryFieldRef()[patchI];
        symmTensorField &alphaPatch = alpha_.boundaryFieldRef()[patchI];
        scalarField &mPatch = m_.boundaryFieldRef()[patchI];
        symmTensorField &FPatch = F_.boundaryFieldRef()[patchI];
        scalarField &KPatch = K_.boundaryFieldRef()[patchI];
        scalarField &GPatch = mu_.boundaryFieldRef()[patchI];
        scalarField &psiPatch = psi_.boundaryFieldRef()[patchI];


        if (debug == 1)
        {
            Info << "calculate boundary stresses" << endl;
        }
        forAll(sigmaP, faceI)
        {
            calculateStress(sigmaP[faceI],
                            dEpsV_Patch[faceI],
                            dEpsD_Patch[faceI],
                            ePatch[faceI],
                            POldPatch[faceI],
                            sOldPatch[faceI],
                            alphaOldPatch[faceI],
                            FOldPatch[faceI],
                            mOldPatch[faceI],
                            dEpsPVPatch[faceI],
                            dEpsPDPatch[faceI],
                            PPatch[faceI],
                            sPatch[faceI],
                            alphaPatch[faceI],
                            mPatch[faceI],
                            FPatch[faceI],
                            KPatch[faceI],
                            GPatch[faceI],
                            psiPatch[faceI],
                            activeYieldP[faceI]);
        }
    }

    // Correct the coupled boundary conditions
    sigma.correctBoundaryConditions();
    activeYield_.correctBoundaryConditions();

    // Under-relax the stress
    sigma.relax();

    if (debug == 1)
    {
        Info << "done calculating new stress" << nl
             << "calculate DSigma_" << endl;
    }
    // Update the effective stress tensor
    // Note: if a_ poro-elasto-plastic law is used then the pore-pressure term
    // will be added after this
    DSigma_ = sigma - sigma.oldTime();
}

void Foam::SANISAND::correct(
    surfaceSymmTensorField &sigma) {}

Foam::scalar
Foam::SANISAND::residual()
{
    // Calculate residual based on change in stress increment
    if (debug == 1)
    {
        Info << "Calculate residual based on change in stress increment" << endl;
    }
    return
        gMax(
            mag(
                DSigma_.primitiveField() - DSigma_.prevIter().primitiveField())) /
        gMax(SMALL + mag(DSigma_.primitiveField()));

}

void Foam::SANISAND::updateTotalFields()
{
    if (debug == 1)
    {
        Info << type() << ": updating total fields" << endl;
    }

    epsilonV_ += DEpsilonV_;
    epsilonD_ += DEpsilonD_;

    DEpsilonP_ = 1.0 / 3.0 * DEpsilonPV_ * I + DEpsilonPD_;
    epsilonP_ += DEpsilonP_;
    r_ = s_ / P_;
    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.primitiveField(), cellI)
    {
        if (activeYield_.primitiveField()[cellI] > SMALL)
        {
            numCellsYielding++;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    Info << "    " << numCellsYielding << " cells are actively yielding"
         << nl << endl;
}



// ************************************************************************* //
