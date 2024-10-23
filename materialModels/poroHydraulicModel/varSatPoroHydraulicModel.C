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

#include "varSatPoroHydraulicModel.H"
#include "volFields.H"


namespace Foam{

    defineTypeNameAndDebug(varSatPoroHydraulicModel, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
varSatPoroHydraulicModel::varSatPoroHydraulicModel(
    const volScalarField &pField,
    const dimensionedVector& gravity)
    : poroHydraulicModel(pField,gravity),
      saturationLawPtr_(pField.mesh().cellZones().size())
{
    forAll(pField.mesh().cellZones(),iZone)
    {
        const dictionary& zoneSubDict =
             poroHydraulicProperties().optionalSubDict
                      (
                        pField.mesh().cellZones().names()[iZone]
                      );

        saturationLawPtr_.set
        (
            iZone,
            saturationLaw::New
            (
                zoneSubDict.get<word>("SWCC"),
                const_cast<dictionary&>(zoneSubDict),
                pField
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

varSatPoroHydraulicModel::~varSatPoroHydraulicModel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> varSatPoroHydraulicModel::pStar() const
{
    tmp<volScalarField> tpStar
    (
        new volScalarField
        (
            IOobject
            (
                "pStar",
                mesh().time().timeName(),
                pField().db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(pField().dimensions(), 0.0),
            "zeroGradient" // Boundaryfield not important for pStar
        )
    );

    volScalarField& pStar_ = tpStar.ref();

    const cellZoneMesh& cellZones = mesh().cellZones();  
    forAll(pStar_,cellI)
    {
        pStar_.internalFieldRef()[cellI] =
            saturationLawPtr_[cellZones.whichZone(cellI)].pStar(cellI);
    }

    pStar_.correctBoundaryConditions();
    
    return tpStar;
}

const tmp<volScalarField> varSatPoroHydraulicModel::kEff(const volScalarField &p)
{
    return k() * kr(p);
}

const tmp<surfaceScalarField> varSatPoroHydraulicModel::kEfff(const volScalarField &p)
{
    return kf()*fvc::interpolate(this->kr(p));
}

tmp<volScalarField> varSatPoroHydraulicModel::kr(const volScalarField &p)
{
    tmp<volScalarField> tkr
    (
        new volScalarField
        (
            IOobject
            (
                "kr",
                mesh().time().timeName(),
                pField().db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );

    volScalarField& kr_ = tkr.ref();

    const cellZoneMesh& cellZones = mesh().cellZones();  
    forAll(p,cellI)
    {
        kr_.internalFieldRef()[cellI] =
            saturationLawPtr_[cellZones.whichZone(cellI)].kr(p.internalField()[cellI],cellI);
    }

    forAll(kr_.boundaryField(),patchI)
    {
        scalarField& krPatch = kr_.boundaryFieldRef()[patchI];
        const scalarField& pPatch = p.boundaryField()[patchI];
        const labelUList& faceCells = mesh().boundaryMesh()[patchI].faceCells();
        forAll(krPatch,faceI)
        {
            label cellI = faceCells[faceI];
            saturationLaw& SModel = saturationLawPtr_[cellZones.whichZone(cellI)];
            krPatch[faceI] = SModel.S(pPatch[faceI],patchI,faceI);
        }
    }
    
    return tkr;
}

tmp<volScalarField> varSatPoroHydraulicModel::S(const volScalarField &p)
{
    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject
            (
                "S",
                mesh().time().timeName(),
                pField().db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );

    volScalarField& S_ = tS.ref();

    const cellZoneMesh& cellZones = mesh().cellZones();  
    forAll(p,cellI)
    {
        S_.internalFieldRef()[cellI] =
            saturationLawPtr_[cellZones.whichZone(cellI)].S(p.internalField()[cellI],cellI);
    }

    forAll(S_.boundaryField(),patchI)
    {
        scalarField& SPatch = S_.boundaryFieldRef()[patchI];
        const scalarField& pPatch = p.boundaryField()[patchI];
        const labelUList& faceCells = mesh().boundaryMesh()[patchI].faceCells();
        forAll(SPatch,faceI)
        {
            label cellI = faceCells[faceI];
            saturationLaw& SModel = saturationLawPtr_[cellZones.whichZone(cellI)];
            SPatch[faceI] = SModel.S(pPatch[faceI],patchI,faceI);
        }
    }
    
    return tS;
}

tmp<volScalarField> varSatPoroHydraulicModel::C(const volScalarField &p)
{
    tmp<volScalarField> tC
    (
        new volScalarField
        (
            IOobject
            (
                "C",
                mesh().time().timeName(),
                pField().db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless/pField().dimensions(), 0.0),
            "zeroGradient" // Boundaryfield not important for C
        )
    );

    volScalarField& C_ = tC.ref();

    const cellZoneMesh& cellZones = mesh().cellZones();  
    forAll(p,cellI)
    {
        C_.internalFieldRef()[cellI] =
            saturationLawPtr_[cellZones.whichZone(cellI)].C(p.internalField()[cellI],cellI);
    }

    C_.correctBoundaryConditions();
    
    return tC;
}



} // End of namespace Foam
// ************************************************************************* //
