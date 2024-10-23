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

#include "poroSolidInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvOption.H"
#include "poroMechanicalLaw2.H"

namespace Foam
{
        defineTypeNameAndDebug(poroSolidInterface, 0);
        defineRunTimeSelectionTable(poroSolidInterface, dictionary);
        addNamedToRunTimeSelectionTable(physicsModel, poroSolidInterface, physicsModel, poroSolid);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::poroSolidInterface::makeIterCtrl()
    {
        iterCtrl_.reset(new iterationControl(const_cast<Time&>(runTime()),interactionProperties_,"Pressure-Displacement"));
    }

void Foam::poroSolidInterface::makePoroFluidCouplingSource()
    {
        dictionary poroSolidSourceDict = subOrEmptyDict("DtoPorePressureCouplingRegions"); //  from v1912 onwards: subDictOrAdd("kinematicPressureSource");
        poroSolidSourceDict.add("type","poroSolidToFluidCouplingSource");
        poroSolidSourceDict.lookupOrAddDefault<word>("selectionMode","all");//  from v1912 onwards: getOrAdd<word>("selectionMode","all");
        poroFluidRef().fvOptions().PtrList<fv::option>::append(fv::option::New("poroSolidToFluidCouplingSource",poroSolidSourceDict,poroFluidMesh()));
    }

bool Foam::poroSolidInterface::checkMechanicalLawUpdateBiotCoeff(const mechanicalLaw& law, const label lawI, PtrList<volScalarField> &bs)
    {
        bool isCoupled = false;
        //Info << law.type() << " " << law.name()  << endl;

        if (law.type()=="poroMechanicalLaw2" || law.type()=="varSatPoroMechanicalLaw")
        {

            const poroMechanicalLaw2 &pmlPtr = dynamic_cast<const poroMechanicalLaw2&>(law);

            bs.set
            (
                lawI,
                new volScalarField(pmlPtr.biotCoeff())
            );

            isCoupled = true;
        }
        else if (law.type()=="poroMechanicalLaw")
        {
            Warning() << " mechanicalLaw 'poroMechanicalLaw' selected,"
                                << " Biot-Willis coefficient will only be regarded in effective stress formulation!" << nl
                                << " Please consider using poroMechanicalLaw2 for consistent treatment!"<<endl;
            isCoupled = true;
            bs.set
            (
                lawI,
                new volScalarField(
                    IOobject(
                    "subMeshBiotCoeff",
                    runTime().timeName(),
                    law.rho()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    ),
                    law.rho()().mesh(),
                    dimensionedScalar("",dimless,1.0)
                    )
            );
        }
        else
        {
            bs.set
            (
                lawI,
                new volScalarField(
                    IOobject(
                    "subMeshBiotCoeff",
                    runTime().timeName(),
                    law.rho()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    ),
                    law.rho()().mesh(),
                    dimensionedScalar("",dimless,1.0)
                    )
            );
            Warning() << "mechanicalLaw " << law.type() << " for material " << law.name()
                            << " is not 'poroMechanicalLaw2' or 'varSatPoroMechanicalLaw'." << nl
                            << "This means there is no coupling from poroFluid to solid, "
                            << "and Biot Coefficient will be set to 1.0 in this region!" << nl
                            << "To account for other values, the storage term in poroHydraulicProperties can be modified."
                            << endl;
        }
        return isCoupled;
    }

void Foam::poroSolidInterface::makeBiotCoeff()
{
    tmp<volScalarField> tb(
        new volScalarField(
            IOobject(
                "biotCoeff",
                runTime().constant(),
                solidMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            solidMesh(),
            dimensionedScalar("",dimless,1.0)
        )
    );
    const PtrList<mechanicalLaw>& laws = dynamic_cast<const PtrList<mechanicalLaw>&>(solid().mechanical());
    // Accumulated subMesh fields and then map to the base mesh
    PtrList<volScalarField> bs(laws.size());

    bool foundPMLaw = false;

    forAll(laws, lawI)
    {
        foundPMLaw = checkMechanicalLawUpdateBiotCoeff(laws[lawI],lawI,bs);
    }
    if(!foundPMLaw)
    {
        WarningInFunction() << "No 'poroMechanicalLaw2' found in 'mechanicalProperties'! " << nl
                            << "This means there is no coupling from poroFluid to solid, "
                            << "and Biot Coefficient will be set to 1.0 everywhere!" << nl
                            << "To account for other values, the storage term in poroHydraulicProperties can be modified."
                            << endl;
    }

    // Map subMesh fields to the base mesh
    if (laws.size()==1)
    {
        tb.ref() = bs[0];
    }
    else
    {
        solid().mechanical().solSubMeshes().mapSubMeshVolFields<scalar>(bs, tb.ref());
    }
    // Clear subMesh fields
    bs.clear();
    if (sharedMesh_)
    {
       b_.reset(tb.ptr());
    }
    else
    {
        b_.reset(solidToPoroFluid().mapTgtToSrc(tb()).ptr());
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poroSolidInterface::poroSolidInterface(
    const word &type,
    Time &runTime,
    const word &region)
    : physicsModel(type, runTime),
      IOdictionary(
          IOobject(
              "poroCouplingProperties",
              runTime.constant(),
              runTime,
              IOobject::MUST_READ,
              IOobject::NO_WRITE)),
      runTime_(runTime),
      interactionProperties_(subDict(type + "Coeffs")),
      sharedMesh_(
          lookupOrAddDefault<Switch>(
              "sharedMesh", true)),
      porosityConstant_(lookupOrAddDefault<Switch>(
              "porosityConstant", true)),
      porosityConstantExplicit_(porosityConstant_
      ? Switch(true)
      : lookupOrAddDefault<Switch>(
              "porosityConstantExplicit", true)),
      solid_(),
      poroFluid_(),
      b_(),
      solidToPoroFluid_(),
      iterCtrl_(),
      intWork_()

{
    solid_ = solidModel::New(runTime, "solid");
    poroFluid_ = poroFluidModel::New(runTime, "poroFluid", sharedMesh_);

    if(!sharedMesh_)
    {
    const word mapMethodName(lookupOrAddDefault<word>("mapMethod","direct"));
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << "unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    Switch consistent(
          lookupOrAddDefault<Switch>(
              "consistent", true));
        if(!consistent)
        {
            if(mapMethodName=="direct")
            {
                WarningInFunction
                    << "direct mapping selected for inconsitent region meshes, changing to imMapNearest"
                    << endl;
                mapMethod = meshToMesh::interpolationMethod::imMapNearest;
            }
        HashTable<word> patchMap;
        wordList cuttingPatches;

        readEntry("patchMap", patchMap);
        readEntry("cuttingPatches", cuttingPatches);
        solidToPoroFluid_.reset
                (
                    new meshToMesh
                    (
                        poroFluidMesh(),
                        solidMesh(),
                        mapMethod,
                        patchMap,
                        cuttingPatches
                    )
                );
        }
        else
        {
        solidToPoroFluid_.reset
                (
                    new meshToMesh
                    (
                        poroFluidMesh(),
                        solidMesh(),
                        mapMethod
                    )
                );
        }
    }

    intWork_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "internalWork",
                runTime.timeName(),
                solidMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh(),
            dimensionedScalar(dimPressure,0.0)
        )
    );

    makePoroFluidCouplingSource();
    Info << "Mesh is shared between Fields: " << sharedMesh_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::poroSolidInterface::~poroSolidInterface()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::poroSolidInterface> Foam::poroSolidInterface::New(
    Time &runTime,
    const word &region)
{
    word fsiTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the poroFluidModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary fsiProperties(
            IOobject(
                "poroCouplingProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));

        fsiProperties.lookup("poroSolidInterface") >> fsiTypeName;
    }

    Info << "Selecting poroSolidInterface " << fsiTypeName << endl;
#if (OPENFOAM >= 2112)
        auto* ctorPtr = dictionaryConstructorTable(fsiTypeName);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                fsiTypeName,
                "poroSolidInterface",
                fsiTypeName,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }
#else
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fsiTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "poroSolidInterface::New(Time&, const word&)")
            << "Unknown poroSolidInterface type " << fsiTypeName
            << endl
            << endl
            << "Valid poroSolidInterface types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

        auto* ctorPtr = cstrIter();
#endif

    return autoPtr<poroSolidInterface>(ctorPtr(runTime, region));
}

void Foam::poroSolidInterface::setDeltaT(Time &runTime)
{
    // For now, the poroFluid sets the time-step
    poroFluidRef().setDeltaT(runTime);
}


void Foam::poroSolidInterface::initializeFields()
{
}


void Foam::poroSolidInterface::movePoroFluidMesh()
{
    notImplemented("additionalPhysics mesh movement not yet implemented for fieldInteractions. Consider shared mesh for moving mesh interactions");
}

void Foam::poroSolidInterface::updateTotalFields()
{}

void Foam::poroSolidInterface::writeFields(const Time &runTime)
{
    autoPtr<volSymmTensorField> DEpsilon;
    if(solid().incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            solidMesh().lookupObject<volTensorField>("grad(DD)");

        DEpsilon.reset(symm(gradDD).ptr());
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            solidMesh().lookupObject<volTensorField>("grad(D)");

        DEpsilon.reset((symm(gradD)-symm(gradD.oldTime())).ptr());
    }

    intWork_.ref() =
        (solid().sigma()-solid().sigma().oldTime())
        && DEpsilon();
    intWork_().write();    

    poroFluidRef().writeFields(runTime);
    solidRef().writeFields(runTime);
}

void Foam::poroSolidInterface::end()
{
    this->IOobject::rename(this->IOobject::name()+".withDefaultValues");
    this->regIOobject::write();
    poroFluidRef().end();
}

// ************************************************************************* //
