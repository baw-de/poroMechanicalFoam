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

#include "anisotropicConstant.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  namespace conductivityModels
  {
    defineTypeNameAndDebug(anisotropicConstant, 0);

    addToRunTimeSelectionTable(
        conductivityModel,
        anisotropicConstant,
        dictionary);

    // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    anisotropicConstant::anisotropicConstant(
        const word &name,
        dictionary &poroHydraulicProperties,
        const volScalarField &pField)
        : conductivityModel(name, poroHydraulicProperties, pField),
          anisotropicConstantProperties(poroHydraulicProperties.subDict(typeName + "Coeffs")),
          anisotropyFactors_(anisotropicConstantProperties.lookup("anisotropyFactors")),
          kf_()
    {
        makeK();

        IOobject anisoHeader(
        "anisotropyFactorField",
        mesh().time().timeName(),
        db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE);


        if(anisoHeader.typeHeaderOk<volSymmTensorField>())
        {
          volSymmTensorField volAnIso(anisoHeader,mesh());
          surfaceScalarField anisotropyFactorField("surfaceNormalAnisotropyFactors",((fvc::interpolate(volAnIso)&mesh().Sf())&mesh().Sf())/pow(mesh().magSf(),2));
          kf_.reset(new surfaceScalarField("kf",fvc::interpolate(k_())*anisotropyFactorField));
          k_.ref() = k_() * 1.0/3.0 * tr(volAnIso);
        }
        else
        {
          surfaceScalarField anisotropyFactorField("surfaceNormalAnisotropyFactors",((anisotropyFactors_&mesh().Sf())&mesh().Sf())/pow(mesh().magSf(),2));
          kf_.reset(new surfaceScalarField("kf",fvc::interpolate(k_())*anisotropyFactorField));
          k_.ref() = k_() * 1.0/3.0 * tr(anisotropyFactors_);
        }
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    tmp<surfaceScalarField> anisotropicConstant::kf() const
    {
      return kf_();
    }

  } // End of namespace conductivityModels
} // End of namespace Foam

//*********************************************************** //
