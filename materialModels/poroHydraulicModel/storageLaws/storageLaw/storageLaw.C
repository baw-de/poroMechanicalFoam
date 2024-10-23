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

#include "storageLaw.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

  defineTypeNameAndDebug(storageLaw, 0);
  defineRunTimeSelectionTable(storageLaw, dictionary);
  // * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

  void storageLaw::makeCs()
  {
    Cs_.reset(new volScalarField(
        IOobject(
          "Cs",
          mesh().time().timeName(),
          db(),
          IOobject::READ_IF_PRESENT, 
          IOobject::NO_WRITE),
        mesh(),
        poroHydraulicProperties_.lookupOrAddDefault<dimensionedScalar>("Cs",dimensionedScalar("", dimless / pField_.dimensions(), 0.0)),
        "zeroGradient"));
  }
  
  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  storageLaw::storageLaw(
      const word &name,
      dictionary &poroHydraulicProperties,
      const volScalarField &pField)
      : name_(name),
        poroHydraulicProperties_(poroHydraulicProperties),
        pField_(pField),
        Cs_(),
        writeSs_(true)
  {}

  // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
