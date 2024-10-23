/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "varSatPoroMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "varSatPoroHydraulicModel.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(varSatPoroMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, varSatPoroMechanicalLaw, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
bool Foam::varSatPoroMechanicalLaw::checkSigmaEffReady(const volSymmTensorField& sigma, const volScalarField& p, const volScalarField& chi)
{
    if (sigmaEff_.valid())
    {
        return true;
    }

    sigmaEff_.reset(
        new volSymmTensorField{
            IOobject
            (
                "sigmaEff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
           sigma + b_ * chi *(p)*symmTensor(I)
        }
    );

    return true;
}

bool Foam::varSatPoroMechanicalLaw::checkSigmaEffReady(const surfaceSymmTensorField& sigma, const surfaceScalarField& p, const surfaceScalarField& chi)
{
    if (sigmaEfff_.valid())
    {
        return true;
    }

    sigmaEfff_.reset(
        new surfaceSymmTensorField{
            IOobject
            (
                "sigmaEfff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
           sigma + b_*chi*(p)*symmTensor(I)
        }
    );
    return true;
}

Foam::effectiveStressModel &Foam::varSatPoroMechanicalLaw::effectiveStressModelRef()
{
    if (!effectiveStressModelPtr_.valid())
    {
        FatalErrorInFunction() << "effectiveStressModel must be initialized first!"
                               << endl;
    }
    return effectiveStressModelPtr_.ref();
}

void Foam::varSatPoroMechanicalLaw::makeN() const
{
    nSubMesh_.reset(
        new const volScalarField(
            IOobject(
                "n_"+mesh().name(),
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dict().get<dimensionedScalar>("n")
        ));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::varSatPoroMechanicalLaw::varSatPoroMechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    poroMechanicalLaw2(name, mesh, dict, nonLinGeom),
    saturationReaderPtr_(readFromDisk_
                         ?new scalarDiskReader("S", mesh, baseMesh(), dict)
                         :NULL),
    effectiveStressModelPtr_(),
    nSubMesh_(),
    rho_water(readFromDisk_
                ? dimensionedScalar(dict.lookupOrDefault<dimensionedScalar>("rho_water",dimensionedScalar("rho_water",dimDensity,1000)))//lookupOrAddDefault<dimensionedScalar>("rho_water",dimensionedScalar("rho_water",dimDensity,1000)))
                : dimensionedScalar("rho_water",dimDensity, 1000)) // Dummy

{
    effectiveStressModelPtr_.reset(effectiveStressModel::New(dict,dict.lookupOrDefault<word>("effectiveStressModel","terzaghi"),mesh));
    if(!buoyancySwitch_)
    {
        Warning << "Buoyancy force is deactivated in variably saturated case!" << nl 
                << "This should only be done, if the weight if the soil is neglected." << nl
                << "For variably saturated calculations please use TOTAL saturated density (bulk density of mixture) and ACTIVATE buoyancy force"
                << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::varSatPoroMechanicalLaw::~varSatPoroMechanicalLaw()
{}

Foam::tmp<Foam::volScalarField> Foam::varSatPoroMechanicalLaw::rho() const
{
    //- Select water density to use
    const dimensionedScalar rhow = readFromDisk_
                                    ? rho_water
                                    : mesh().time().lookupObject<varSatPoroHydraulicModel>("poroHydraulicModel").rho();


    //- Check if saturation should be read from file (this would be higherst priority if true)
    const volScalarField SRef = readFromDisk_
                                ? saturationReader().field()
                                : lookupFluidField("S");

    const volScalarField nRef = readFromDisk_
                                ? n()
                                : lookupFluidField("n");

    tmp<volScalarField> rho_total(
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            rhoScalar()-(1-SRef)*nRef*rhow,
            "zeroGradient"
        )
    );

    return rho_total;
}

void Foam::varSatPoroMechanicalLaw::correct(volSymmTensorField& sigma)
{

    //- Check if saturation should be read from file (this would be higherst priority if true)
    tmp<volScalarField> SRef = readFromDisk_
                                ? saturationReader().field()
                                : lookupFluidField("S");


    tmp<volScalarField> nRef = readFromDisk_
                                ? n()
                                : lookupFluidField("n");

    // Lookup the pressure field
    tmp<volScalarField> pRef = readFromDisk_
                                ? pReader().field()
                                : lookupFluidField(pName_);

    tmp<volScalarField> chiRef(effectiveStressModelPtr_->chi(nRef(),SRef(),pRef()));

    // check if sigmaEff has been initialized (should be done only once per calculation)
    this->checkSigmaEffReady(sigma,pRef(),chiRef());

    // Calculate effective stress
    //-- Note that we could just pass "sigma" here but we use a separate field
    //-- called sigmaEff just for post-processing visualisation of the effective
    //-- stress <-- for stress state dependent material laws like Mohr-Coulomb it is important to use sigmaEff since the strength depends on tr(sigmaEff)
    effectiveStressMechLawPtr_->correct(sigmaEff_());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure

    sigma = sigmaEff_() - b_*chiRef()*(pRef())*symmTensor(I);
}


void Foam::varSatPoroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    //- Check if saturation should be read from file (this would be higherst priority if true)
    tmp<volScalarField> SRef = readFromDisk_
                                ? saturationReader().field()
                                : lookupFluidField("S");


    tmp<volScalarField> nRef = readFromDisk_
                                ? n()
                                : lookupFluidField("n");

    // Lookup the pressure field
    tmp<volScalarField> pRef = readFromDisk_
                                ? pReader().field()
                                : lookupFluidField(pName_);

    // Interpolate pressure to the faces
    const surfaceScalarField pf(fvc::interpolate(pRef));

    // Interpolate pressure to the faces
    const surfaceScalarField Sf(fvc::interpolate(SRef));

    const surfaceScalarField chif = fvc::interpolate(effectiveStressModelPtr_->chi(nRef,SRef,pRef));

    // check if sigmaEff has been initialized (should be done only once per calculation)
    checkSigmaEffReady(sigma,pf,chif);

    // Calculate effective stress
    effectiveStressMechLawPtr_->correct(sigmaEfff_());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma = sigmaEfff_() - b_*chif*(pf)*symmTensor(I);
}


// ************************************************************************* //
