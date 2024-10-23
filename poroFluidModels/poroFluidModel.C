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

#include "poroFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(poroFluidModel, 0);
defineRunTimeSelectionTable(poroFluidModel, dictionary);
addNamedToRunTimeSelectionTable(physicsModel, poroFluidModel, physicsModel,
                                poroFluid);

// * * * * * * * * * * * * * * Privat Member Functions * * * * * * * * * * * * *
// //
void poroFluidModel::makeIterCtrl()
{
    iterCtrl_.reset
    (
        new iterationControl
        (
            const_cast<Time&>(runTime()),
            poroFluidDict(),
             "Pressure"
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

poroFluidModel::poroFluidModel
(
    const word& type,
    Time& runTime,
    const word& fieldName,
    const word& region,
    const bool sharedMesh
)
    : physicsModel(type, runTime),
      dynamicFvMesh
      (
        IOobject
        (
            region,
            bool(region == dynamicFvMesh::defaultRegion)
                ? fileName(runTime.caseConstant())
                : sharedMesh
                    ? fileName(runTime.caseConstant() / region)
                    : fileName(runTime.caseConstant() / "solid"),
            runTime, IOobject::MUST_READ, IOobject::NO_WRITE
        )
       ),
      name_(type),
      sharedMesh_(sharedMesh), // Is 'no' by default but can be set
                               // to yes for coupled calculations
      runTime_(runTime),
      poroFluidProperties_
      (
        IOobject
        (
            "poroFluidProperties",
            runTime.constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
      ),
      type_(type),
      fieldName_(fieldName),
      g_
      (
        IOobject // will be initialized in the constructor body
        (
            "g",
            runTime.constant(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector::zero)
      ),
      nPtr_(),
      n0Ptr_(),
      p_
      (
        IOobject
        (
            "p",
            runTime.timeName(),
                  // if mesh is shared between solid and addPhysics register
                  // fields in this objectRegistery instead of (solid-) mesh
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("", dimPressure, 0.0)
      ),
      pDot_
      (
        IOobject
        (
            "pDot",
            runTime.timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(p_)
      ),
      pRGHheader_("p_rgh", runTime.timeName(), *this, IOobject::MUST_READ),
      p_rgh_
      (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
                      // if mesh is shared between solid and addPhysics register
                      // fields in this objectRegistery instead of (solid-) mesh
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("", dimPressure, 0.0)
      ),
      hydraulicGradient_
      (
        IOobject
        (
            "i",
            runTime.timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("", dimless, vector::zero)
      ),
      phi_
      (
        IOobject
        (
            "phi",
            runTime.timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimVolume / dimTime, 0.0)
      ),
      fvOptions_((fv::options::New(mesh()))),
      coeffsDictName_(type_ + "Coeffs"),
      iterCtrl_(),
      maxDeltaT_
      (
        runTime.controlDict().lookupOrDefault<scalar>
        (
                       "maxDeltaT", GREAT
        )
      ),
      CoNumber_(runTime.controlDict().lookupOrDefault<scalar>("CoNumber", 1))
//               relaxationMethod_(
//                   poroFluidDict().lookupOrDefault<Foam::word>("relaxationMethod",
//                   "fixed"))
{

    // Check grad(p) into objectRegistry. Why does this field not checkIn
    // automatically? (Denis)
    hydraulicGradient_.checkIn();
    mesh().setFluxRequired(fieldName_);

    IOobject gheader("g", runTime.constant(), *this, IOobject::MUST_READ);
    if (!runTime.objectRegistry::foundObject<fvMesh>("solid")
        || gheader.typeHeaderOk<uniformDimensionedVectorField>())
    {
        Info << "reading gravity for poroFluid from file" << endl;
        g_ = uniformDimensionedVectorField(gheader);
        if (mag(g_).value() == 0.0)
        {
            FatalErrorIn("poroFluidModel::poroFluidModel(const "
                         "word&,Time&,const word&,const word&,const bool)")
                << "Gravity is off, this leads to problems in the poroFluid "
                   "calculation."
                << nl
                << "Please provide a non-zero-valued gravity file for "
                   "poroFluid!"
                << exit(FatalError);
        }
    }
    else
    {
        Info << "reading solid gravity for poroFluid" << endl;
        g_ = runTime.subRegistry("solid")
                 .objectRegistry::lookupObject<uniformDimensionedVectorField>(
                     "g");
        if (mag(g_).value() == 0.0)
        {
            FatalErrorIn("poroFluidModel::poroFluidModel(const "
                         "word&,Time&,const word&,const word&,const bool)")
                << "Solid gravity is off, this leads to problems in the "
                   "poroFluid calculation."
                << nl
                << "In case you want to keep solid gravity off (e.g. for "
                   "weightless calculations),"
                << "please provide a seperate non-zero-valued gravity file for "
                   "poroFluid!"
                << exit(FatalError);
        }
    }


    Info << "poroFluidModel: " << type << nl
         << "with solution field: " << fieldName_ << nl
         << "on mesh: " << p_.mesh().name() << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<poroFluidModel>
poroFluidModel::New(Time& runTime, const word& region, const bool sharedMesh)
{
    word poroFluidModelTypeName;

    // In case a case that is set up for a coupled calculation is supposed to be
    // run as simple (uncoupled) poroFluid calculation (i. e. to check the pure
    // poroFluid convergence and set up) the solver can be set to poroFluid and
    // in the controlDict under the subDict 'poroFluid' one can change the
    // region to 'poroFluid'. This also allows the region name of a coupled
    // calulation not to be 'poroFluid' but user defined.
    word runRegion(runTime.controlDict()
                       .subOrEmptyDict("poroFluid")
                       .lookupOrDefault<word>("region", region));

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid model is created, otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary poroFluidProperties(
            IOobject("poroFluidProperties",
                     runRegion == dynamicFvMesh::defaultRegion
                         ? fileName(runTime.caseConstant())
                         : fileName(runTime.caseConstant() + "/" + runRegion),
                     runTime, IOobject::MUST_READ, IOobject::NO_WRITE));
        poroFluidProperties.lookup("poroFluidModel") >> poroFluidModelTypeName;
    }

    //- for compatebility
    if (poroFluidModelTypeName == "extPoroFluid")
    {
        poroFluidModelTypeName = "varSatPoroFluid";
    }

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(poroFluidModelTypeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup(poroFluidModelTypeName, "poroFluidModel",
                             poroFluidModelTypeName,
                             *dictionaryConstructorTablePtr_)
            << exit(FatalIOError);
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(poroFluidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("poroFluidModel::New(Time&, const fvMesh&, const bool) ")
            << "Unknown poroFluidModel type " << poroFluidModelTypeName << endl
            << endl
            << "Valid  poroFluidModels are : " << endl
            << dictionaryConstructorTablePtr_->toc() << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif
    return autoPtr<poroFluidModel>(ctorPtr(runTime, runRegion, sharedMesh));
}

void Foam::poroFluidModel::pRGHisRequired()
{
    if (!pRGHheader_.typeHeaderOk<volScalarField>(true))
    {
        FatalErrorIn(type() + "::pRGHisRequired()")
            << "This poroFluidModel requires the 'p_rgh' field to be specified!"
            << abort(FatalError);
    }
}

const volScalarField& Foam::poroFluidModel::n()
{
    //- If not yet happend, initialize n with n0
    if(!nPtr_.valid())
    {
        nPtr_.reset(n0().ptr());
    }
    //Otherwise we return the field
    return nPtr_();
}

void Foam::poroFluidModel::update_porosity(const volScalarField Dn, const bool incremental)
{
    // porosity is changing, and we dont have a field for the changing porosity yet
    // we initalize with n0 (the porosity inside poroHydraulicModel)
    if(!n0Ptr_.valid())
    {
        n0Ptr_.reset
        (
            new volScalarField
            (
                "n0",
                n()
            )
        );   
    }

    if(incremental)
    {
        nPtr_.ref() = n().oldTime() + Dn;
    }
    else
    {
        nPtr_.ref() = n0Ptr_() + Dn;
    }
}

void Foam::poroFluidModel::writeFields(const Time& runTime)
{
    mesh().objectRegistry::write();
}

bool Foam::poroFluidModel::update()
{
    notImplemented("moving poroFluidMesh is not yet implemented");
    return false;
}

void Foam::poroFluidModel::end()
{
    poroFluidProperties_.IOobject::rename(poroFluidProperties_.IOobject::name()
                                          + ".withDefaultValues");
    poroFluidProperties_.regIOobject::write();
    physicsModel::end();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
