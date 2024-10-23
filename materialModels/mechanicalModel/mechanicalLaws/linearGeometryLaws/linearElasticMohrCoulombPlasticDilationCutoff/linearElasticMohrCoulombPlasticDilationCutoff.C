/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-exsigmaBtend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: yousigmaB can redistribute it and/or modify it
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

#include "linearElasticMohrCoulombPlasticDilationCutoff.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMohrCoulombPlasticDilationCutoff, 0);
    addToRunTimeSelectionTable(mechanicalLaw, linearElasticMohrCoulombPlasticDilationCutoff, linGeomMechLaw);
}

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

void Foam::linearElasticMohrCoulombPlasticDilationCutoff::calculateStress(
    symmTensor &sigma,
    scalar &activeYield,
    scalar &f,
    scalar &ePore,
    scalar &EigenError) const
{
    // Principal stresses
    vector sigma_prin = vector::zero;

    // Principal directions
    tensor ev = tensor::zero;

    activeYield = 0;

    for (int i = 0; i < 6; i++)
    {
        if (sigma[i] > SMALL)
        {
            sigma[i] += (i + 1) * (1e-6) * sigma[i];
        }
    }

    calculateEigens(sigma_prin, ev, sigma, EigenError);

    // Re-order the principal stresses
    scalar sigma1 = min(sigma_prin[0], 0);
    scalar sigma2 = min(sigma_prin[1], 0);
    scalar sigma3 = min(sigma_prin[2], 0);

    label sigma1_po = 0;
    label sigma2_po = 1;
    label sigma3_po = 2;

    for (int i = 1; i < 3; i++)
    {
        if (sigma_prin[i] > sigma1)
        {
            sigma1 = min(sigma_prin[i], 0);
            sigma1_po = i;
        }
    }

    for (int i = 1; i >= 0; i--)
    {
        if (sigma_prin[i] < sigma3)
        {
            sigma3 = min(sigma_prin[i], 0);
            sigma3_po = i;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if ((i != sigma1_po) && (i != sigma3_po))
        {
            sigma2_po = i;
            sigma2 = min(sigma_prin[i], 0);
        }
    }

    vector sigmaB = vector(sigma1, sigma2, sigma3);

    vector b_(b0_);
    vector rp_(rp0_);
    vector r_lg1_(r_lg10_);
    vector r_lg2_(r_lg20_);

    // Evaluate the yielding function f
    f = k_ * sigmaB[0] - sigmaB[2] - 2 * c_.value() * Foam::sqrt(k_);

    // Check if yielding
    if (f > SMALL)
    {

        activeYield = 1.0;
        // Determine the return type
        if (ePore >= eCrit_.value())
        {
            activeYield = 2.0;

            b_ = vector(1, 0, -1);
            rp_ = vector((C_ & b_) / (b_ & (C_ & a_)));
            r_lg1_ = vector(1, 1, 1);
            r_lg2_ = vector(1, 1, 1);
        }
        // Calculate the boundary planes
        const scalar t1 =
            (r_lg1_ & (invC_ & (sigma_prin - sigma_a_))) / (r_lg1_ & (invC_ & r_lf1_));

        const scalar t2 =
            (r_lg2_ & (invC_ & (sigma_prin - sigma_a_))) / (r_lg2_ & (invC_ & r_lf2_));

        const scalar p_12 = (rp_ ^ r_lf1_) & (sigma_prin - sigma_a_);

        const scalar p_13 = (rp_ ^ r_lf2_) & (sigma_prin - sigma_a_);

        // Compute the stress field based on the correct return type

        if ((p_12 >= 0) && (p_13 <= 0))
        {
            sigmaB -= f * rp_;
        }
        else if ((p_12 < 0) && (p_13 < 0))
        {
            sigmaB = t1 * r_lf1_ + sigma_a_;
        }
        else if ((p_12 > 0) && (p_13 > 0))
        {
            sigmaB = t2 * r_lf2_ + sigma_a_;
        }
        else if ((t1 > 0) && (t2 > 0))
        {
            sigmaB = sigma_a_;
        }

        sigma_prin[sigma1_po] = min(sigmaB[0], 0);
        sigma_prin[sigma2_po] = min(sigmaB[1], 0);
        sigma_prin[sigma3_po] = min(sigmaB[2], 0);

        // Form the diagonal stress
        const diagTensor sigma_diag =
            diagTensor(sigma_prin[0], sigma_prin[1], sigma_prin[2]);

        // Transform the returned stress back into general space
        sigma = symm(ev.T() & sigma_diag & ev);
    }
}

void Foam::linearElasticMohrCoulombPlasticDilationCutoff::calculateEigens(
    vector &sigma_prin,
    tensor &ev,
    const symmTensor sigma,
    scalar &EigenError) const
{
    // Check for a zero stress tensor
    if (mag(sigma) < SMALL)
    {
        ev = I;

        return;
    }

    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if (
        (
            mag(sigma.xy()) + mag(sigma.xz()) + mag(sigma.xy()) + mag(sigma.yz()) + mag(sigma.xz()) + mag(sigma.yz())) < 1e-6 * (mag(sigma.xx()) + mag(sigma.yy()) + mag(sigma.zz())))
    {
        // Diagonal matrix
        i = sigma.xx();
        ii = sigma.yy();
        iii = sigma.zz();
    }
    else
    {
        const scalar a = -sigma.xx() - sigma.yy() - sigma.zz();

        const scalar b =
            sigma.xx() * sigma.yy() + sigma.xx() * sigma.zz() + sigma.yy() * sigma.zz() - sigma.xy() * sigma.xy() - sigma.xz() * sigma.xz() - sigma.yz() * sigma.yz();

        const scalar c =
            -sigma.xx() * sigma.yy() * sigma.zz() - sigma.xy() * sigma.yz() * sigma.xz() - sigma.xz() * sigma.xy() * sigma.yz() + sigma.xz() * sigma.yy() * sigma.xz() + sigma.xy() * sigma.xy() * sigma.zz() + sigma.xx() * sigma.yz() * sigma.yz();

        // If there is a zero root
        if (mag(c) < SMALL)
        {
            const scalar disc = Foam::max(sqr(a) - 4 * b, 0.0);
            if (debug)
            {
                WarningIn("poroMohrCoulob::calculateEigens(...)")
                    << "Stress tensor has a zero root! symTensor: " << sigma << endl;
            }
            const scalar q = -0.5 * Foam::sqrt(max(scalar(0), disc));

            i = 0;
            ii = -0.5 * a + q;
            iii = -0.5 * a - q;

            EigenError = 1;
        }
        else
        {
            const scalar Q = (a * a - 3.0 * b) / 9.0;
            const scalar R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

            const scalar R2 = sqr(R);
            const scalar Q3 = pow3(Q);

            // Three different real roots
            if (R2 < Q3)
            {
                const scalar sqrtQ = Foam::sqrt(Q);
                const scalar theta = Foam::acos(R / (Q * sqrtQ));

                const scalar m2SqrtQ = -2 * sqrtQ;
                const scalar aBy3 = a / 3;

                i = m2SqrtQ * Foam::cos(theta / 3) - aBy3;

                ii =
                    m2SqrtQ * Foam::cos((theta + constant::mathematical::twoPi) / 3.0) - aBy3;
                iii =
                    m2SqrtQ * Foam::cos((theta - constant::mathematical::twoPi) / 3.0) - aBy3;

            }
            else
            {
                const scalar A = Foam::cbrt(R + Foam::sqrt(R2 - Q3));

                // Three equal real roots
                if (A < SMALL)
                {
                    const scalar root = -a / 3;
                    i = root;
                    ii = root;
                    iii = root;
                }
                else
                {
                    if (debug)
                    {
                        // Complex roots
                        WarningIn("poroMohrCoulob::calculateEigens(...)")
                            << "Complex eigenvalues detected for symmTensor: "
                            << sigma << nl << "Setting roots to zero!" << endl;
                        i = 0;
                        ii = 0;
                        iii = 0;
                    }
                    EigenError = 2;
                }
            }
        }
    }

    // Sort the eigenvalues into ascending order
    if (mag(sigma.xy()) < 1e-10 && mag(sigma.xz()) < 1e-10 && mag(sigma.yz()) < 1e-10)
    {
        ev = I;
        i = sigma.xx();
        ii = sigma.yy();
        iii = sigma.zz();
        if (debug)
        {
            WarningIn("poroMohrCoulob::calculateEigens(...)")
                << "almost-eigenform detected for symmTensor: "
                << sigma << nl << "Setting roots to diagonal terms!" << endl;
        }
        EigenError = 3;
    }

    if (mag(i) > mag(ii))
    {
        Swap(i, ii);
    }

    if (mag(ii) > mag(iii))
    {
        Swap(ii, iii);
    }

    if (mag(i) > mag(ii))
    {
        Swap(i, ii);
    }

    sigma_prin[0] = i;
    sigma_prin[1] = ii;
    sigma_prin[2] = iii;

    if (mag(sigma.xy()) < 1e-10 && mag(sigma.xz()) < 1e-10 && mag(sigma.yz()) < 1e-10)
    {
        if (debug)
        {
            Info << "returning" << endl;
        }
        return;
    }

    for (int j = 0; j < 3; j++)
    {
        if (mag(sigma_prin[j]) < SMALL)
        {
            if (j == 0)
            {
                ev[j * 3 + 0] = 1;
                ev[j * 3 + 1] = 0;
                ev[j * 3 + 2] = 0;
            }
            else if (j == 1)
            {
                ev[j * 3 + 0] = 0;
                ev[j * 3 + 1] = 1;
                ev[j * 3 + 2] = 0;
            }
            else
            {
                ev[j * 3 + 0] = 0;
                ev[j * 3 + 1] = 0;
                ev[j * 3 + 2] = 1;
            }
        }
        else
        {
            const symmTensor A = symmTensor(sigma - sigma_prin[j] * I);

            // Calculate the sub-determinants of the 3 components
            const scalar sd0 = A.yy() * A.zz() - A.yz() * A.yz();
            const scalar sd1 = A.xx() * A.zz() - A.xz() * A.xz();
            const scalar sd2 = A.xx() * A.yy() - A.xy() * A.xy();

            const scalar magSd0 = mag(sd0);
            const scalar magSd1 = mag(sd1);
            const scalar magSd2 = mag(sd2);

            // Evaluate the eigenvector using the largest sub-determinant
            if ((magSd0 > magSd1) && (magSd0 > magSd2) && (magSd0 > SMALL))
            {
                vector newEv =
                    vector(
                        1,
                        (A.yz() * A.xz() - A.zz() * A.xy()) / sd0,
                        (A.yz() * A.xy() - A.yy() * A.xz()) / sd0);
                newEv /= mag(newEv);

                ev[j * 3] = newEv[0];
                ev[j * 3 + 1] = newEv[1];
                ev[j * 3 + 2] = newEv[2];
            }
            else if ((magSd1 > magSd2) && (magSd1 > SMALL))
            {
                vector newEv =
                    vector(
                        (A.xz() * A.yz() - A.zz() * A.xy()) / sd1,
                        1,
                        (A.xz() * A.xy() - A.xx() * A.yz()) / sd1);
                newEv /= mag(newEv);

                ev[j * 3] = newEv[0];
                ev[j * 3 + 1] = newEv[1];
                ev[j * 3 + 2] = newEv[2];
            }
            else if (magSd2 > SMALL)
            {
                vector newEv =
                    vector(
                        (A.xy() * A.yz() - A.yy() * A.xz()) / sd2,
                        (A.xy() * A.xz() - A.xx() * A.yz()) / sd2,
                        1);
                newEv /= mag(newEv);

                ev[j * 3] = newEv[0];
                ev[j * 3 + 1] = newEv[1];
                ev[j * 3 + 2] = newEv[2];
            }
            else
            {
                if (mag(sigma.xy()) < SMALL && mag(sigma.xz()) < SMALL && mag(sigma.yz()) < SMALL)
                {
                    ev = I;
                    i = sigma.xx();
                    ii = sigma.yy();
                    iii = sigma.zz();
                    if (debug)
                    {
                        WarningIn("poroMohrCoulob::calculateEigens(...)")
                            << "almost-eigenform detected for symmTensor: "
                            << sigma << nl << "Setting roots to diagonal terms!" << endl;
                    }
                    EigenError = 3;
                }
                else
                {
                    if (debug)
                    {
                        WarningIn("poroMohrCoulob::calculateEigens(...)")
                            << "Strange things detected for stress tensor: "
                            << sigma << nl
                            << "Setting eigenvectors to (0, 0, 0)" << endl;
                    }
                    ev[j * 3] = 0;
                    ev[j * 3 + 1] = 0;
                    ev[j * 3 + 2] = 0;
                    EigenError = 4;
                }
            }
        }
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticMohrCoulombPlasticDilationCutoff::linearElasticMohrCoulombPlasticDilationCutoff(
    const word &name,
    const fvMesh &mesh,
    const dictionary &dict,
    const nonLinearGeometry::nonLinearType &nonLinGeom)
    : mechanicalLaw(name, mesh, dict, nonLinGeom),
      //rho_(dict.lookup("rho")),
      E_(dict.get<dimensionedScalar>("E")),
      nu_(dict.get<dimensionedScalar>("nu")),
      lambda_(
          planeStress()
              ? nu_ * E_ / ((1.0 + nu_) * (1.0 - nu_))
              : nu_ * E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_))),
      mu_(E_ / (2.0 * (1.0 + nu_))),
      K_(lambda_ + (2.0 / 3.0) * mu_),
      varPhi_(dict.get<dimensionedScalar>("frictionAngle")),
      c_(dict.get<dimensionedScalar>("cohesion")),
      varPsi0_(dict.get<dimensionedScalar>("dilationAngle")),
      e0_(dict.get<dimensionedScalar>("e0")),
      eCrit_(dict.get<dimensionedScalar>("e_crit")),
      e_(
          IOobject(
              "e",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          e0_),
      k_(
          (
              (1 + sin(varPhi_ / 180.0 * constant::mathematical::pi)) / (1 - sin(varPhi_ / 180.0 * constant::mathematical::pi)))
              .value()
              ),
      m0_(
          (
              (1 + sin(varPsi0_ / 180.0 * constant::mathematical::pi)) / (1 - sin(varPsi0_ / 180.0 * constant::mathematical::pi)))
              .value()
              ),
      a_(k_, 0, -1),
      b0_(m0_, 0, -1),
      C_(
          (K_ + (4.0 / 3.0) * mu_).value(),
          (K_ - (2.0 / 3.0) * mu_).value(),
          (K_ - (2.0 / 3.0) * mu_).value(),
          (K_ + (4.0 / 3.0) * mu_).value(),
          (K_ - (2.0 / 3.0) * mu_).value(),
          (K_ + (4.0 / 3.0) * mu_).value()),
      invC_(inv(C_)),
      rp0_((C_ & b0_) / (b0_ & (C_ & a_))),
      r_lf1_(1, 1, k_),
      r_lf2_(1, k_, k_),
      r_lg10_(1, 1, m0_),
      r_lg20_(1, m0_, m0_),
      sigma_a_((2.0 * c_ * Foam::sqrt(k_) / (k_ - 1.0)).value() * vector::one),
      Dsigma_(
          IOobject(
              "Dsigma",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)),
      Dsigmaf_(
          IOobject(
              "Dsigmaf",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)),
      DEpsilonP_(
          IOobject(
              "DEpsilonP",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedSymmTensor("zero", dimless, symmTensor::zero)),
      DEpsilonPf_(
          IOobject(
              "DEpsilonPf",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimensionedSymmTensor("zero", dimless, symmTensor::zero)),
      epsilonP_(
          IOobject(
              "epsilonP",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedSymmTensor("zero", dimless, symmTensor::zero)),
      epsilonPEq_(
          IOobject(
              "epsilonPEq",
              mesh.time().timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("zero", dimless, 0.0)),
      activeYield_(
          IOobject(
              "activeYield",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("0", dimless, 0)),
      f_(
          IOobject(
              "yieldFunction",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("0", dimless, 0)),
      EigenError_(
          IOobject(
              "EigenError",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("0", dimless, 0))
{
    e_.storeOldTime();
    Info << "Mohr-Coulomb DilationCap Parameters:" << nl
         << "Cohesion: " << c_ << nl
         << "Friction Angle: " << varPhi_ << nl
         << "Dilation Angle: " << varPsi0_ << nl
         << "Critical Pore Ratio: " << eCrit_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticMohrCoulombPlasticDilationCutoff::~linearElasticMohrCoulombPlasticDilationCutoff()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Now located in mechanicalLaw
/*Foam::tmp<Foam::volScalarField>
Foam::linearElasticMohrCoulombPlasticDilationCutoff::rho() const
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
}*/

Foam::tmp<Foam::volScalarField>
Foam::linearElasticMohrCoulombPlasticDilationCutoff::impK() const
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
            2.0 * mu_ + lambda_));
}

const Foam::dimensionedScalar &Foam::linearElasticMohrCoulombPlasticDilationCutoff::mu() const
{
    return mu_;
}

const Foam::dimensionedScalar &
Foam::linearElasticMohrCoulombPlasticDilationCutoff::lambda() const
{
    return lambda_;
}

void Foam::linearElasticMohrCoulombPlasticDilationCutoff::correct(volSymmTensorField &sigma)
{
    if(gSum(mag(sigma0())())==0)
    {
        sigma0()=sigma.oldTime();
    }

    // Update epsilon
    updateEpsilon();

    // Update the increment of strain
    volSymmTensorField DEpsilon(epsilon() - epsilon().oldTime());

    // Calculate trial elastic stress assuming Hooke's elastic law
    const volSymmTensorField DSigmaTrial
    (
        2.0*mu_*DEpsilon + lambda_*I*tr(DEpsilon)
    );

    // update pore ratio
    e_ = e_.oldTime() + tr(DEpsilon);

    // Set sigma to sigma effective trial, including the initial stress
    sigma = Dsigma_.oldTime() + DSigmaTrial + sigma0();

    // Take a reference to internal fields for efficiency
    symmTensorField &sigmaI = sigma.primitiveFieldRef();
    scalarField &fI = f_.primitiveFieldRef();
    scalarField &eI = e_.primitiveFieldRef();
    scalarField &EigenErrorI = EigenError_.primitiveFieldRef();
    scalarField &activeYieldI = activeYield_.primitiveFieldRef();


    // Correct sigma internal field
    forAll(sigmaI, cellI)
    {
        calculateStress(sigmaI[cellI], activeYieldI[cellI], fI[cellI], eI[cellI], EigenErrorI[cellI]);
    }

    // Correct sigma on the boundary patches
    forAll(sigma.boundaryField(), patchI)
    {
        // Take references to the boundary patches for efficiency
        symmTensorField &sigmaP = sigma.boundaryFieldRef()[patchI];
        scalarField &fP = f_.boundaryFieldRef()[patchI];
        scalarField &eP = e_.boundaryFieldRef()[patchI];
        scalarField &EigenErrorP = EigenError_.boundaryFieldRef()[patchI];
        scalarField &activeYieldP = activeYield_.boundaryFieldRef()[patchI];


        if (!sigma.boundaryField()[patchI].coupled())
        {
            forAll(sigmaP, faceI)
            {
                calculateStress(sigmaP[faceI], activeYieldP[faceI], fP[faceI], eP[faceI], EigenErrorP[faceI]);
            }
        }
    }

    // Correct the coupled boundary conditions
    sigma.correctBoundaryConditions();
    activeYield_.correctBoundaryConditions();

    // Store previous iteration of Dsigma as it is used to calculate the
    // residual
    Dsigma_.storePrevIter();

    // Update the stress variation tensor
    Dsigma_ = sigma - sigma0();
    // Note: if a poro-elasto-plastic law is used then the pore-pressure term
    // will be added to the effective stress after this
}
void Foam::linearElasticMohrCoulombPlasticDilationCutoff::correct
(
    surfaceSymmTensorField& sigma
)
{
    NotImplemented;
    /*
    // Update epsilon
    updateEpsilonf();

    // Update the increment of strain
    DEpsilonf_ = epsilonf() - epsilonf().oldTime();

    // Calculate trial elastic stress assuming Hooke's elastic law
    const surfaceSymmTensorField DSigmaTrial
    (
        2.0*mu_*DEpsilonf_ + lambda_*I*tr(DEpsilonf_)
    );

    // update pore ratio
    e_ = e_.oldTime() + tr(DEpsilon);

    // Set sigma to sigma effective trial, including the initial stress
    sigma = Dsigma_.oldTime() + DSigmaTrial + sigma0();

    // Take a reference to internal fields for efficiency
    symmTensorField& sigmaI = sigma;
    scalarField& activeYieldI = activeYield_;

#ifdef OPENFOAM_COM
    const labelList& faceOwner = mesh().faceOwner();
    const labelList& faceNeighbour = mesh().faceNeighbour();
#else
    const unallocLabelList& faceOwner = mesh().faceOwner();
    const unallocLabelList& faceNeighbour = mesh().faceNeighbour();
#endif

    // Correct sigma internal field
    forAll(sigmaI, faceI)
    {
        const label ownCellID = faceOwner[faceI];
        const label neiCellID = faceNeighbour[faceI];

        calculateStress(sigmaI[faceI], activeYieldI[ownCellID]);

        // Update the neighbour activeYield as it is a vol field, i.e. if a face
        // is yielding then we set the active yield flag for both the owner and
        // neighbour cells
        activeYieldI[neiCellID] = activeYieldI[ownCellID];
    }

    // Correct sigma on the boundary patches
    forAll(sigma.boundaryField(), patchI)
    {
        // Take references to the boundary patches for efficiency
#ifdef OPENFOAM_NOT_EXTEND
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryFieldRef()[patchI];
#else
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryField()[patchI];
#endif

        forAll(sigmaP, faceI)
        {
            calculateStress(sigmaP[faceI], activeYieldP[faceI]);
        }
    }

    // Store previous iteration of Dsigma as it is used to calculate the
    // residual
    Dsigmaf_.storePrevIter();

    // Update the stress variation tensor
    Dsigmaf_ = sigma - sigma0f();
    // Note: if a poro-elasto-plastic law is used then the pore-pressure term
    // will be added after this
    */
}

Foam::scalar Foam::linearElasticMohrCoulombPlasticDilationCutoff::residual()
{
    // Calculate residual based on change in plastic strain increment
    if (
        mesh().foundObject<surfaceVectorField>("grad(D)f") || mesh().foundObject<surfaceVectorField>("grad(DD)f"))
    {
        return
            gMax(
                mag(
                    Dsigmaf_.primitiveField() - Dsigmaf_.prevIter().primitiveField())) /
            gMax(SMALL + mag(Dsigmaf_.primitiveField()));
    }
    /*     else
    {
        return
#ifdef OPENFOAMESIORFOUNDATION
            gMax(
                mag(
                    Dsigma_.primitiveField() - Dsigma_.prevIter().primitiveField())) /
            gMax(SMALL + mag(Dsigma_.primitiveField()));
#else
            gMax(
                mag(
                    Dsigma_.internalField() - Dsigma_.prevIter().internalField())) /
            gMax(SMALL + mag(Dsigma_.internalField()));
#endif
    } */
    else
    {
        scalarField x = (mag(Dsigma_.primitiveField() - Dsigma_.prevIter().primitiveField())) /
                        (SMALL + (mag(sigma0().primitiveField())));
        return pow(
            gSum(
                pow(x,
                    2)),
            0.5);
    }
}

/**
 * Updates the total fields of the linearElasticMohrCoulombPlasticDilationCutoff
 * class. This function calculates the increment of plastic strain and the
 * increment of total strain. It also updates the total plastic strain and the
 * total equivalent plastic strain. Additionally, it counts the number of
 * cells actively yielding.
 *
 * @throws None
 */
void Foam::linearElasticMohrCoulombPlasticDilationCutoff::updateTotalFields()
{
    Info << type() << ": updating total fields" << endl;

    // Calculate increment of plastic strain
    // This is complicated because DEpsilonP also has a volumetric term
    symmTensorField &DEpsilonPI = DEpsilonP_.primitiveFieldRef();

    // Calculate the increment of total strain
    volSymmTensorField DEpsilon(
        IOobject(
            "DEpsilon",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero));

    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField &gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        DEpsilon = symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField &gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate gradient of displacement increment
        const volTensorField gradDD = gradD - gradD.oldTime();

        DEpsilon = symm(gradDD);
    }

    const symmTensorField &DEpsilonI = DEpsilon.primitiveField();
    const symmTensorField &sigmaEffI = Dsigma_.primitiveField();

    const scalar mu = mu_.value();
    const scalar lambda = lambda_.value();

    forAll(DEpsilonPI, cellI)
    {
        // Off-diagonal strains
        DEpsilonPI[cellI].xy() =
            DEpsilonI[cellI].xy() - (sigmaEffI[cellI].xy() / (2.0 * mu));
        DEpsilonPI[cellI].xz() =
            DEpsilonI[cellI].xz() - (sigmaEffI[cellI].xz() / (2.0 * mu));
        DEpsilonPI[cellI].yz() =
            DEpsilonI[cellI].yz() - (sigmaEffI[cellI].yz() / (2.0 * mu));

        // Solve a linear system (Ax = b) to calculate the strains on the
        // diagonal strains
        const tensor A =
            tensor(
                -(2.0 * mu + lambda), -lambda, -lambda,
                -lambda, -(2.0 * mu + lambda), -lambda,
                -lambda, -lambda, -(2.0 * mu + lambda));

        const vector b =
            vector(
                sigmaEffI[cellI].xx() - (2.0 * mu + lambda) * DEpsilonI[cellI].xx() - lambda * DEpsilonI[cellI].yy() - lambda * DEpsilonI[cellI].zz(),

                sigmaEffI[cellI].yy() - lambda * DEpsilonI[cellI].xx() - (2.0 * mu + lambda) * DEpsilonI[cellI].yy() - lambda * DEpsilonI[cellI].zz(),

                sigmaEffI[cellI].zz() - lambda * DEpsilonI[cellI].xx() - lambda * DEpsilonI[cellI].yy() - (2.0 * mu + lambda) * DEpsilonI[cellI].zz());

        const vector x = (inv(A) & b);

        DEpsilonPI[cellI].xx() = x.x();
        DEpsilonPI[cellI].yy() = x.y();
        DEpsilonPI[cellI].zz() = x.z();
    }

    forAll(DEpsilonP_.boundaryField(), patchI)
    {
        symmTensorField &DEpsilonPP = DEpsilonP_.boundaryFieldRef()[patchI];

        const symmTensorField &DEpsilonP = DEpsilon.boundaryField()[patchI];
        const symmTensorField &sigmaEffP = Dsigma_.boundaryField()[patchI];

        forAll(DEpsilonPP, faceI)
        {
            // Off-diagonal strains
            DEpsilonPP[faceI].xy() =
                DEpsilonP[faceI].xy() - (sigmaEffP[faceI].xy() / (2.0 * mu));
            DEpsilonPP[faceI].xz() =
                DEpsilonP[faceI].xz() - (sigmaEffP[faceI].xz() / (2.0 * mu));
            DEpsilonPP[faceI].yz() =
                DEpsilonP[faceI].yz() - (sigmaEffP[faceI].yz() / (2.0 * mu));

            // Solve a linear system (Ax = b) to calculate the strains on the
            // diagonal strains
            const tensor A =
                tensor(
                    -(2.0 * mu + lambda), -lambda, -lambda,
                    -lambda, -(2.0 * mu + lambda), -lambda,
                    -lambda, -lambda, -(2.0 * mu + lambda));

            const vector b =
                vector(
                    sigmaEffP[faceI].xx() - (2.0 * mu + lambda) * DEpsilonP[faceI].xx() - lambda * DEpsilonP[faceI].yy() - lambda * DEpsilonP[faceI].zz(),

                    sigmaEffP[faceI].yy() - lambda * DEpsilonP[faceI].xx() - (2.0 * mu + lambda) * DEpsilonP[faceI].yy() - lambda * DEpsilonP[faceI].zz(),

                    sigmaEffP[faceI].zz() - lambda * DEpsilonP[faceI].xx() - lambda * DEpsilonP[faceI].yy() - (2.0 * mu + lambda) * DEpsilonP[faceI].zz());

            const vector x = (inv(A) & b);

            DEpsilonPP[faceI].xx() = x.x();
            DEpsilonPP[faceI].yy() = x.y();
            DEpsilonPP[faceI].zz() = x.z();
        }
    }

    //     DEpsilon
    //   - (1.0/3.0)*I*tr(sigma - Dsigma_.oldTime())/(2.0*mu_ + 3.0*lambda_)
    //   - dev(sigma - Dsigma_.oldTime())/(2.0*mu_);

    // Calculate the increment of equivalent plastic strain
    const volScalarField DEpsilonPEq = sqrt((2.0 / 3.0) * magSqr(dev(DEpsilonP_)));

    // Update the total plastic strain
    epsilonP_ = (epsilonP_.oldTime() + DEpsilonP_);

    // Update the total equivalent plastic strain
    epsilonPEq_ = (epsilonPEq_.oldTime() + DEpsilonPEq);

    Info << "    Max DEpsilonPEq is " << gMax(DEpsilonPEq) << endl;
    Info << "    Max epsilonPEq is " << gMax(epsilonPEq_) << endl;

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
