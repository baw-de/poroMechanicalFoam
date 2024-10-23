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

#include "MFront.H"
#include <MGIS/Behaviour/Integrate.hxx>
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MFront, 0);
    addToRunTimeSelectionTable(
        mechanicalLaw, MFront, linGeomMechLaw);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*const Foam::word *Foam::MFront::toString(mgis::behaviour::Behaviour::Kinematic kin)
{
    switch (kin)
    {
    case mgis::behaviour::Behaviour::Kinematic::UNDEFINEDKINEMATIC:
        return "UNDEFINEDKINEMATIC";
    case mgis::behaviour::Behaviour::Kinematic::SMALLSTRAINKINEMATIC:
        return "SMALLSTRAINKINEMATIC";
    case mgis::behaviour::Behaviour::Kinematic::COHESIVEZONEKINEMATIC:
        return "COHESIVEZONEKINEMATIC";
    case mgis::behaviour::Behaviour::Kinematic::FINITESTRAINKINEMATIC_F_CAUCHY:
        return "FINITESTRAINKINEMATIC_F_CAUCHY";
    case mgis::behaviour::Behaviour::Kinematic::FINITESTRAINKINEMATIC_ETO_PK1:
        return "FINITESTRAINKINEMATIC_ETO_PK1";
    }

    FatalErrorIn("Foam::MFront::toString(mgis::behaviour::Behaviour::Kinematic kin) unknown kinematic");
}
const Foam::word *Foam::MFront::toString(mgis::behaviour::Behaviour::Symmetry sym)
{
    switch (sym)
    {
    case mgis::behaviour::Behaviour::Symmetry::ISOTROPIC:
        return "ISOTROPIC";
    case mgis::behaviour::Behaviour::Symmetry::ORTHOTROPIC:
        return "ORTHOTROPIC";
    }
    FatalErrorIn("Foam::MFront::toString(mgis::behaviour::Behaviour::Symmetry sym) unknown symmetry");
}
const Foam::word *Foam::MFront::btypeToString(int btype)
{
    if (btype == mgis::behaviour::Behaviour::GENERALBEHAVIOUR)
        return "GENERALBEHAVIOUR";
    if (btype == mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR)
        return "STANDARDSTRAINBASEDBEHAVIOUR";
    if (btype == mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR)
        return "STANDARDFINITESTRAINBEHAVIOUR";
    if (btype == mgis::behaviour::Behaviour::COHESIVEZONEMODEL)
        return "COHESIVEZONEMODEL";

    FatalErrorIn("Foam::MFront::btypeToString(int btype)");
}
const Foam::word *Foam::MFront::varTypeToString(int v)
{
    if (v == mgis::behaviour::Variable::SCALAR)
        return "SCALAR";
    if (v == mgis::behaviour::Variable::VECTOR)
        return "VECTOR";
    if (v == mgis::behaviour::Variable::STENSOR)
        return "STENSOR";
    if (v == mgis::behaviour::Variable::TENSOR)
        return "TENSOR";
    FatalErrorIn("Foam::MFront::varTypeToString(int v).");
}*/

// Construct from dictionary
Foam::MFront::MFront(
    const word &name,
    const fvMesh &mesh,
    const dictionary &dict,
    const nonLinearGeometry::nonLinearType &nonLinGeom)
    : mechanicalLaw(name, mesh, dict, nonLinGeom),
      MFrontLibPath(dict.lookup("MFrontLibPath")),
      behaviourName(dict.lookup("MFrontBehaviour")),
      hypothesis_(mgis::behaviour::Hypothesis::TRIDIMENSIONAL),
      behaviour_(mgis::behaviour::load(MFrontLibPath, behaviourName, hypothesis_)),
      properties_(dict.subDict("MFrontProperties")),
      rho_(dict.lookup("rho")),
      epsilon_(
          IOobject(
              "epsilon",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimensionedSymmTensor("zero", dimless, symmTensor::zero)),
      MFrontBehaviourData(behaviour_, int(epsilon_.internalField().size()))
{
    // Force storage of strain old time
    epsilon_.oldTime();
    if (behaviour_.gradients.size() != 1)
        FatalErrorIn(
            "The behaviour must have exactly a single gradient as input.");

    if (behaviour_.gradients[0].name != "Strain")
        FatalErrorIn("The behaviour must be driven by strain.");

    if (behaviour_.gradients[0].type != mgis::behaviour::Variable::STENSOR)
        FatalErrorIn("Strain must be a symmetric tensor.");

    /*if (mgis::behaviour::getVariableSize(_behaviour.gradients[0], hypothesis) !=
         MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime)
         FatalErrorIn("Strain must have {:d} components instead of {:d}.",
                   MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime,
                   mgis::behaviour::getVariableSize(_behaviour.gradients[0],
                                                    hypothesis));*/

    if (behaviour_.thermodynamic_forces.size() != 1)
        FatalErrorIn(
            "The behaviour must compute exactly one thermodynamic force.");

    if (behaviour_.thermodynamic_forces[0].name != "Stress")
        FatalErrorIn("The behaviour must compute stress.");

    if (behaviour_.thermodynamic_forces[0].type !=
        mgis::behaviour::Variable::STENSOR)
        FatalErrorIn("Stress must be a symmetric tensor.");

    /* if (mgis::behaviour::getVariableSize(_behaviour.thermodynamic_forces[0],
                                          hypothesis) !=
         MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime)
         FatalErrorIn("Stress must have {:d} components instead of {:d}.",
                   MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime,
                   mgis::behaviour::getVariableSize(
                       _behaviour.thermodynamic_forces[0], hypothesis));*/

    if (!behaviour_.esvs.empty())
    {
        notImplemented("No external state variable support implemented");
        /*if (_behaviour.esvs[0].name != "Temperature")
         {
             FatalErrorIn(
                 "Only temperature is supported as external state variable.");
         }*/
    }

    /*if (behaviour_.mps.size() != properties_.size())
    {
        FatalErrorIn("Number of material properties in MFront does not match number of passed elements");
    }*/
    Info << "Behaviour:    " << behaviour_.behaviour << endl;
    Info << "Hypothesis:   " << mgis::behaviour::toString(hypothesis_) << endl;
    Info << "Source:       " << behaviour_.source << endl;
    Info << "TFEL version: " << behaviour_.tfel_version << endl;
    //INFO("Behaviour type: `{:s}'.", btypeToString(behaviour.btype));
    //INFO("Kinematic:      `{:s}'.", toString(behaviour.kinematic));
    //INFO("Symmetry:       `{:s}'.", toString(behaviour.symmetry));

    updateParams();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MFront::~MFront()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MFront::updateParams()
{
    if (!behaviour_.mps.empty())
    {
        for (uint iParam = 0; iParam < behaviour_.mps.size(); iParam++)
        {
            if (behaviour_.mps[iParam].type == mgis::behaviour::Variable::STENSOR || behaviour_.mps[iParam].type == mgis::behaviour::Variable::TENSOR)
            {
                notImplemented("Foam::MFront::updateParams() Tensorial material properties are not yet implemented");
            }
            scalar tmpParam(readScalar(properties_.lookup(behaviour_.mps[iParam].name)));
            Info << "Setting MFront Parameter" << behaviour_.mps[iParam].name << " = " << tmpParam << endl;
            mgis::behaviour::setParameter(behaviour_, behaviour_.mps[iParam].name, tmpParam);
        }
    }
}

/*Foam::symmTensor Foam::MFront::MFrontToOf(mgis::behaviour::Variable::SYMMTENSOR &sT){
    return ()}*/

Foam::tmp<Foam::volScalarField> Foam::MFront::rho() const
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

#ifdef OPENFOAMESIORFOUNDATION
    tresult.ref().correctBoundaryConditions();
#else
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}

Foam::tmp<Foam::volScalarField> Foam::MFront::bulkModulus() const
{
    notImplemented("Foam::MFront::bulkModulus()");

    // Kepp the compiler happy
    return rho();
}

const Foam::dimensionedScalar &Foam::MFront::rhoScalar() const
{
    return rho_;
}

Foam::tmp<Foam::volScalarField> Foam::MFront::impK() const
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
            dimensionedScalar("", dimPressure, 1.0)));
}

const Foam::dimensionedScalar &Foam::MFront::mu() const
{
    notImplemented("Foam::MFront::mu()");

    // Kepp the compiler happy
    return rhoScalar();
}

const Foam::dimensionedScalar &Foam::MFront::K() const
{
    notImplemented("Foam::MFront::K()");

    // Kepp the compiler happy
    return rhoScalar();
}

const Foam::dimensionedScalar &Foam::MFront::E() const
{
    notImplemented("Foam::MFront::E()");

    // Kepp the compiler happy
    return rhoScalar();
}

const Foam::dimensionedScalar &Foam::MFront::nu() const
{
    notImplemented("Foam::MFront::nu()");

    // Kepp the compiler happy
    return rhoScalar();
}

const Foam::dimensionedScalar &Foam::MFront::lambda() const
{
    notImplemented("Foam::MFront::lambda()");

    // Kepp the compiler happy
    return rhoScalar();
}

void Foam::MFront::correct(volSymmTensorField &sigma)
{
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField &gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon_ = epsilon_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField &gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon_ = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn(
                "void Foam::MFront::"
                "correct(volSymmTensorField& sigma)")
                << "Not implemented for planeStress!" << abort(FatalError);
        }
    }
    //auto const &eps_m_prev =
    symmTensorField &epsI = epsilon_.internalField();
    forAll(epsI, iCell)
    {
        symmTensor &epsITemp = epsI[iCell];
        std::vector<double> epsIArray = {epsITemp.XX, epsITemp.YY, epsITemp.ZZ, epsITemp.XY, epsITemp.XZ, epsITemp.YZ};
        std::copy_n(epsIArray.data(), 6,
                    MFrontBehaviourData.s0.gradients.data() + MFrontBehaviourData.s0.gradients_stride);
    }

    mgis::behaviour::integrate(MFrontBehaviourData, behaviour_)
}

void Foam::MFront::correct(surfaceSymmTensorField &sigma)
{
    notImplemented("Foam::MFront::correct(surfaceSymmTensorField)");
}

// ************************************************************************* //
