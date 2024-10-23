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

#include "poroMechanicalLaw2.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroMechanicalLaw2, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, poroMechanicalLaw2, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
bool Foam::poroMechanicalLaw2::checkSigmaEffReady(const volSymmTensorField& sigma, const volScalarField& p)
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
           sigma + b_*(p)*symmTensor(I)
        }
    );
    return true;
}

bool Foam::poroMechanicalLaw2::checkSigmaEffReady(const surfaceSymmTensorField& sigma, const surfaceScalarField& p)
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
           sigma + b_*(p)*symmTensor(I)
        }
    );
    return true;
}

Foam::tmp<Foam::volScalarField> Foam::poroMechanicalLaw2::lookupFluidField(const word fieldName) const
{
    tmp<volScalarField> tScalar;
    if(mesh()==baseMesh())
    {
        tScalar.reset(
            new volScalarField
                (
                    fieldName + "_copy",
                    baseMesh().lookupObject<volScalarField>
                        (
                            fieldName
                        )
                )
            );
    }
    else
    {
        tScalar.reset(
             baseMesh().lookupObject<mechanicalModel>
            (
            "mechanicalProperties"
            ).solSubMeshes().lookupBaseMeshVolField<scalar>
            (fieldName,mesh()));
    }
    return tScalar;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::poroMechanicalLaw2::poroMechanicalLaw2
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    effectiveStressMechLawPtr_
    (
        mechanicalLaw::NewLinGeomMechLaw
        (
            dict.subDict("effectiveStressMechanicalLaw").get<word>("type"),
            mesh,
            dict.subDict("effectiveStressMechanicalLaw"),
            nonLinGeom
        )
    ),
    sigmaEff_(),
    sigmaEfff_(),
    buoyancySwitch_(
        dict.lookupOrDefault<Switch>//lookupOrAddDefault<Switch>
        (
            "buoyancy", true
        )
    ),
    pName_(dict.lookupOrDefault<word>("pressureFieldName", buoyancySwitch_?"p":"p_rgh")),//lookupOrAddDefault<word>("pressureFieldName", buoyancySwitch_?"p":"p_rgh")),
    readFromDisk_(dict.lookupOrDefault<Switch>("readPressureFromDisk", false)),//lookupOrAddDefault<Switch>("readPressureFromDisk", false)),
    pReaderPtr_(readFromDisk_
                ?new scalarDiskReader(pName_, mesh, baseMesh(), dict)
                :NULL
    ),
    b_
    (
        dict.lookupOrDefault<dimensionedScalar>//lookupOrAddDefault<dimensionedScalar>
        (
            "biotCoeff", dimensionedScalar("0", dimless, 1.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poroMechanicalLaw2::~poroMechanicalLaw2()
{}


Foam::tmp<Foam::volScalarField> Foam::poroMechanicalLaw2::impK() const
{
    return effectiveStressMechLawPtr_->impK();
}


void Foam::poroMechanicalLaw2::correct(volSymmTensorField& sigma)
{
    // Lookup the pressure field
    tmp<volScalarField> p = readFromDisk_
                                ? pReader().field()
                                : lookupFluidField(pName_);

    // check if sigmaEff has been initialized (should be done only once per calculation)
    checkSigmaEffReady(sigma,p());

    // Calculate effective stress
    //-- Note that we could just pass "sigma" here but we use a separate field
    //-- called sigmaEff just for post-processing visualisation of the effective
    //-- stress <-- for stress state dependent material laws like Mohr-Coulomb it is important to use sigmaEff since the strength depends on tr(sigmaEff)
    effectiveStressMechLawPtr_->correct(sigmaEff_());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma == sigmaEff_() - b_*(p())*symmTensor(I);
}


void Foam::poroMechanicalLaw2::correct(surfaceSymmTensorField& sigma)
{
    // Lookup the pressure field
    tmp<volScalarField> pRef = readFromDisk_
                                ? pReader().field()
                                : lookupFluidField(pName_);

    // Interpolate pressure to the faces
    const surfaceScalarField pf(fvc::interpolate(pRef()));

    // check if sigmaEff has been initialized (should be done only once per calculation)
    checkSigmaEffReady(sigma,pf);

    // Calculate effective stress
    effectiveStressMechLawPtr_->correct(sigmaEfff_());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma = sigmaEfff_() - b_*(pf)*symmTensor(I);
}


// ************************************************************************* //
