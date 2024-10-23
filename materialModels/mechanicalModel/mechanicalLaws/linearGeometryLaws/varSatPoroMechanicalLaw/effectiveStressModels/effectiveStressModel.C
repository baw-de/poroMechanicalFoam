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

#include "effectiveStressModel.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  defineTypeNameAndDebug(effectiveStressModel, 0);
  defineRunTimeSelectionTable(effectiveStressModel, dictionary);
  // * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

  void effectiveStressModel::makeChi()
  {
    chi_.reset(
                new volScalarField(
                IOobject(
                    "chi",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("",dimless,1.0)
            )
      );
  }

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  effectiveStressModel::effectiveStressModel(
      const dictionary &effectiveStressModelDict,
      const word &name,
      const fvMesh &mesh)
      : name_(name),
        effectiveStressModelDict_(effectiveStressModelDict),
        chi_(),
        mesh_(mesh)
  {
  }

  tmp<volScalarField> effectiveStressModel::chi(const volScalarField &n, const volScalarField &S, const volScalarField &p)
  {

    if(!chi_.valid())
    {
      makeChi();
    }

    // Calculate internal field
    chi_.ref().primitiveFieldRef() = chi(n.internalField(), S.internalField(), p.internalField())();

    // Calculate boundary fields
    forAll(chi_().boundaryField(),iPatch)
    {
      if(chi_().boundaryField().types()[iPatch]!="empty")
      {
        chi_.ref().boundaryFieldRef()[iPatch] = chi(n.boundaryField()[iPatch], S.boundaryField()[iPatch], p.boundaryField()[iPatch])();
      }
    }

    return chi_();
  }
  // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
