/*---------------------------------------------------------------------------*\
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

#include "ohdeElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "pointFieldFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ohdeElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, ohdeElastic, linGeomMechLaw
    );
}


void Foam::ohdeElastic::updateStiffness(const volSymmTensorField& sigma)
{

    scalar oneThird = 1./3.;

    Info << "Calculating current K, mu" << endl;
    K_ == K0_* pow(oneThird*tr(sigma)/p0_,n_);
    mu_ == mu0_* pow(oneThird*tr(sigma)/p0_,n_);

    Info << "Current Kmax:" << Foam::max(K_) << nl
	 << "Current Kmin:" << Foam::min(K_) << nl
	 << "Current mu_max:" << Foam::max(mu_) << nl
	 << "Current mu_min:" << Foam::min(mu_) << nl
	 << "Calculating current E, nu" << endl;
    // Calculate Young's modulus
    E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

    // Calculate Poisson's ratio
    nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));

    Info << "Calculating current lamda" << endl;
    // Set first Lame parameter
    if (planeStress())
        {
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));

            if (solvePressureEqn())
            {
                FatalErrorIn
                (
                    "ohdeElasticMisesPlastic::ohdeElasticMisesPlastic::()"
                )   << "planeStress must be 'off' when solvePressureEqn is "
                    << "enabled" << abort(FatalError);
            }
        }
        else
        {
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::ohdeElastic::ohdeElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu0_("mu", dimPressure, 0.0),
    K0_("K0", dimPressure, 0.0),
    n_("n", dimless, 0.0),
    p0_("p0", dimPressure, 0.0),
    K_
    (
        IOobject
	    (
		"K" + name,
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("0", dimPressure, 0.0)
    ),
    mu_
    (
        IOobject
	    (
		"mu_" + name,
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("0", dimPressure, 0.0)
    ),
    E_
    (
        IOobject
	    (
		"E_" + name,
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("0", dimPressure, 0.0)
    ),
    nu_
    (
        IOobject
	    (
		"nu_" + name,
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("0", dimless, 0.0)
    ),
    lambda_
    (
        IOobject
	    (
		"lambda_" + name,
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("0", dimPressure, 0.0)
    )
{
    // Store old times
    epsilon().storeOldTime();
    epsilonf().storeOldTime();

    // Read shear modulus
    mu0_ = dict.get<dimensionedScalar>("mu0");

    // Read bulk modulus
    K0_ = dict.get<dimensionedScalar>("K0");

    // Read shear modulus
    p0_ = dict.get<dimensionedScalar>("p0");

    // Read bulk modulus
    n_ = dimensionedScalar(dict.lookup("n"));

    // Read the initial stress
    if (dict.found("sigma0"))
    {
        Info<< "Reading sigma0 from the dict" << endl;
        sigma0() = dimensionedSymmTensor(dict.lookup("sigma0"));
    }
    else if (gMax(mag(sigma0())()) > SMALL)
    {
        Info<< "Reading sigma0 stress field" << endl;
    }

    updateStiffness(sigma0());

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::ohdeElastic::bulkModulus() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            K_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().correctBoundaryConditions();
#else
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::ohdeElastic::impK() const
{
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "impK",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                2.0*mu_ + lambda_
            )
        );
}


void Foam::ohdeElastic::correct(volSymmTensorField& sigma)
{
    // Update epsilon
    updateEpsilon();

    // Calculate stress using epsilon
    correct(sigma, epsilon());

    updateStiffness(sigma + sigma0());
}


void Foam::ohdeElastic::correct
(
    volSymmTensorField& sigma,
    const volSymmTensorField& epsilon
)
{
    // Hooke's law: partitioned deviatoric and dilation form
    sigma = 2.0*mu_*dev(epsilon) + K_*tr(epsilon)*I ;
}


// ************************************************************************* //
