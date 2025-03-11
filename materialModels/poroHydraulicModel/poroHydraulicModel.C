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

#include "poroHydraulicModel.H"
#include "volFields.H"


namespace Foam{
    defineTypeNameAndDebug(poroHydraulicModel, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
poroHydraulicModel::poroHydraulicModel
(
    const volScalarField& pField,
    const dimensionedVector& gravity
)
    : regIOobject
    (
        IOobject
        (
            "poroHydraulicModel",
            pField.time().constant(),
            pField.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    poroHydraulicProperties_
    (
        IOobject
        (
            "poroHydraulicProperties",
            pField.time().constant(),
            pField.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    pField_(pField),
    storageLawPtr_(pField.mesh().cellZones().size()),
    conductivityModelPtr_(pField.mesh().cellZones().size()),
    rho_
    (
        poroHydraulicProperties_.lookupOrAddDefault<dimensionedScalar>
        (
            "rho",
            dimensionedScalar("rho", dimDensity, 1000)
        )
    ),
    gamma_
    (
        IOobject
        (
            "gamma_water",
            mesh().time().constant(),
            pField.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*gravity
    ),
    href_
    (
        poroHydraulicProperties_.lookupOrAddDefault<dimensionedScalar>
        (
            "href",
            dimensionedScalar("href", dimLength, 0.0)
        )
    ),
    z_
    ( //  geodetic height (depends on direction of gravity)
        IOobject
        (
            "z",
            mesh().time().constant(),
            pField.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        - mesh().C() & vector(gamma_.value()).normalise()
    ),
    p_Hyd_
    ( //  hydrostatic pressure
        IOobject
        (
            "p_Hyd",
            mesh().time().timeName(),
            pField.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (href_ - z_) * magGamma()
    ),
    updatesSs_(false),
    updatesK_(false)
{

    Info << "Creating the poroHydraulicModel" << nl
         << "Gravity direction: " << vector(gamma_.value()).normalise() << nl
         << "Water specific weight: " << magGamma() << nl
	 << "Referential Watertable is at z = " << href_ << endl;

    forAll(pField.mesh().cellZones(),zoneI)
    {
        const dictionary& zoneSubDict =
             poroHydraulicProperties_.optionalSubDict 
                      (
                        pField.mesh().cellZones().names()[zoneI]
                      );

        storageLawPtr_.set
        (
            zoneI,
            storageLaw::New
            (
                zoneSubDict.get<word>("StorageModel"),
                const_cast<dictionary&>(zoneSubDict),
                pField_
            )
        );

        updatesSs_ = (updatesSs_ || storageLawPtr_[zoneI].updatesSs());

        conductivityModelPtr_.set
        (
            zoneI,
            conductivityModel::New
            (
                zoneSubDict.getOrDefault<word>("conductivityModel","constant"),
                const_cast<dictionary&>(zoneSubDict),
                pField
            )
        );

        updatesK_ = (updatesK_ || conductivityModelPtr_[zoneI].updatesK());
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

poroHydraulicModel::~poroHydraulicModel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> poroHydraulicModel::n0() const
{
    IOobject nHeader
    (
        "n",
        mesh().time().constant(),
        pField().db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (nHeader.typeHeaderOk<volScalarField>(true))
    {
        tmp<volScalarField> tn0
        (
            new volScalarField
            (
                nHeader,
                mesh()
            )
        );

        return tn0;
    }
    else
    {
        tmp<volScalarField> tn0
        (
            new volScalarField
            (
                IOobject
                (
                    "n",
                    mesh().time().constant(),
                    pField().db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar(dimless,0.0)
            )
        );

        volScalarField& n = tn0.ref();
        
        const cellZoneMesh& cellZones = mesh().cellZones();
        PtrList<dimensionedScalar> nList(cellZones.size());  

	if(cellZones.size()==0)
	{
		FatalError << "No cellZones detected, please add a cellZone to your case!"
			   << endl;
	}

        forAll(cellZones,iZone)
        {
            const dictionary& zoneSubDict =
                    poroHydraulicProperties().optionalSubDict 
                    (
                        cellZones.names()[iZone]
                    );
            nList.set
            (
                iZone,
                new dimensionedScalar(zoneSubDict.get<dimensionedScalar>("n"))
            );
        }

        forAll(n,cellI)
        {
            n.internalFieldRef()[cellI] = nList[cellZones.whichZone(cellI)].value();
        }

        forAll(n.boundaryField(),patchI)
        {
            scalarField& nPatch = n.boundaryFieldRef()[patchI];
            const labelUList& faceCells = mesh().boundaryMesh()[patchI].faceCells();
            forAll(nPatch,faceI)
            {
                label cellI = faceCells[faceI];
                nPatch[faceI] = nList[cellZones.whichZone(cellI)].value();
            }
        }

        return tn0;
    }
}

tmp<volScalarField> poroHydraulicModel::Ss(const volScalarField& n, const volScalarField& p)
{
    tmp<volScalarField> tSs
    (
        new volScalarField
        (
            IOobject
            (
                "Ss",
                mesh().time().timeName(),
                pField().db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless / pField().dimensions() , 0.0),
            "zeroGradient" // Boundaryfield not important for Ss
        )
    );

    volScalarField& Ss_ = tSs.ref();

    const cellZoneMesh& cellZones = mesh().cellZones();  
    forAll(p,cellI)
    {
        Ss_.internalFieldRef()[cellI] =
            storageLawPtr_[cellZones.whichZone(cellI)].Ss
                (
                    n.internalField()[cellI],
                    p.internalField()[cellI],
                    cellI
                );
    }

    Ss_.correctBoundaryConditions();

    return tSs;
}

const tmp<volScalarField> poroHydraulicModel::k() const
{
    tmp<volScalarField> tk
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                mesh().time().timeName(),
                pField().db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimLength / dimTime, 0.0)
        )
    );

    volScalarField& k_ = tk.ref();

    const cellZoneMesh& cellZones = mesh().cellZones();  
    forAll(k_,cellI)
    {
        conductivityModel& cModel = conductivityModelPtr_[cellZones.whichZone(cellI)];
        k_.internalFieldRef()[cellI] = cModel.k(cellI);
    }

    forAll(k_.boundaryField(),patchI)
    {
        scalarField& kPatch = k_.boundaryFieldRef()[patchI];
        const labelUList& faceCells = mesh().boundaryMesh()[patchI].faceCells();
        forAll(kPatch,faceI)
        {
            label cellI = faceCells[faceI];
            conductivityModel& cModel = conductivityModelPtr_[cellZones.whichZone(cellI)];
            kPatch[faceI] = cModel.k(patchI,faceI);
        }
    }
    
    return tk;
}

const tmp<surfaceScalarField> poroHydraulicModel::kf() const
{
    tmp<volScalarField> tk = this->k();
    return fvc::interpolate(tk());
}

} // End of namespace Foam
// ************************************************************************* //
