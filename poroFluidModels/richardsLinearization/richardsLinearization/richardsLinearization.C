/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "richardsLinearization.H"
#include "fvc.H"
#include "fvm.H"
#include "fvMatrix.H"
//#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  defineTypeNameAndDebug(richardsLinearization, 0);
  defineRunTimeSelectionTable(richardsLinearization, dictionary);

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  richardsLinearization::richardsLinearization(
      const word &name,
      varSatPoroHydraulicModel &poroHydraulic,
      dictionary &poroFluidProperties,
      volScalarField& S)
      : name_(name),
        poroHydraulic_(poroHydraulic),
        S_(S),
        C_
        (
          IOobject
          (
            "C",
            S.time().timeName(),
            S.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          S_.mesh(),
          dimensionedScalar(dimless/poroHydraulic.pField().dimensions(),0.0)
        ),
        mesh_(poroHydraulic_.mesh())
  {
  }

  // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

  scalar richardsLinearization::L2Norm(const scalarField &x) const
  {
    return pow(
        gSum(
            pow(x, 2)),
        0.5);
  }

  /*scalar richardsLinearization::MassBalance(const volScalarField &pField, const volScalarField &n, const volScalarField &S, const volScalarField &flux) const
  {
    volScalarField theta("theta", n * S);
    scalarField volMassChangeField = volMassChange(theta);
    scalarField miscStorageVolume = volMiscStorage(pField, n);
    scalarField volFlux = flux * mesh().V();

    if (pField.time().objectRegistry::foundObject<volScalarField>(word("nDot")))
    {
      Info << "Coupled Storage detected" << endl;
      const volScalarField &nDot = pField.time().objectRegistry::lookupObject<volScalarField>("nDot");
      const volScalarField &chi = pField.time().objectRegistry::lookupObject<volScalarField>("chi");
      #ifdef OPENFOAMESIORFOUNDATION
      miscStorageVolume += (chi * nDot).ref().primitiveField();
      #else
      miscStorageVolume += (chi * nDot)().internalField();
      #endif
    }

    volScalarField MassBalance_(
        IOobject(
            "MassBalanceError",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh(),
        dimensionedScalar("",dimless,0.0)
    );
    #ifdef OPENFOAMESIORFOUNDATION
    MassBalance_.primitiveFieldRef() = mag(volMassChangeField + miscStorageVolume - volFlux);
    #else
    MassBalance_.internalField() = mag(volMassChangeField + miscStorageVolume - volFlux);
    #endif
    // Info << "Volume: " << sum(mesh().V()) << endl;

    scalar volMassChange_ = gSum(volMassChangeField);
    scalar miscStorage = gSum(miscStorageVolume);
    scalar fluxBoundary = gSum(volFlux);

    scalar massBalanceError =
        (volMassChange_ + miscStorage - fluxBoundary);

    // Detailed Info too expensive for now, deactivated (infoFreqency might be a solution)
    /*Info << "| "
         << "Water Volume Change" << tab << tab << "| " << volMassChange_ << tab << "m³/s" << endl
         << "| "
         << "Miscellaneous Storage" << tab << tab << "| " << miscStorage << tab << "m³/s" << endl
         << "| "
         << "Boundary Flux" << tab << tab << tab << "| " << fluxBoundary << tab << "m³/s" << endl
         << "| "
         << "rel. Mass Balance Error" << tab << tab << "| " << massBalanceError << tab << "%" << endl
         << "-------------------------------------------------------------------" << endl;
    */
  /*  return massBalanceError;
  }

  tmp<scalarField> richardsLinearization::volMassChange(const volScalarField &theta) const
  {
    tmp<scalarField> Vw(
        new scalarField(fvc::ddt(theta) * mesh().V()));
    return Vw;
  }

  tmp<scalarField> richardsLinearization::volMiscStorage(const volScalarField &pField, const volScalarField &n) const
  {
    tmp<scalarField> MStor(
        new scalarField(pField.size(), 0.0));
    if (name() != "steadyState")
    {
    #ifdef OPENFOAMESIORFOUNDATION
      scalarField &MStor_ = MStor.ref();
    #else
      scalarField &MStor_ = MStor();
    #endif

      const volScalarField &Ss = poroHydraulic().Ss();
      MStor_ = (n * fvc::ddt(Ss, pField)) * mesh().V();
    }
    return MStor;
  }

  tmp<volScalarField> richardsLinearization::massFlux(const surfaceScalarField &k_eff, const volScalarField &pField) const
  {
    tmp<volScalarField> massFlux(
        new volScalarField(
            IOobject(
                "massFlux",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            (fvc::laplacian(k_eff, pField))));
    return massFlux;
  }*/

  /*tmp<fvScalarMatrix> richardsLinearization::ddpk(const surfaceScalarField &kField, volScalarField &pField)
  {
    tmp<fvScalarMatrix> tkMatrix(
      new fvScalarMatrix(
        pField,
        kField.dimensions() * pField.dimensions() * dimVol / dimArea));
        return tkMatrix;
  }*/

} // namespace Foam

// ************************************************************************* //
