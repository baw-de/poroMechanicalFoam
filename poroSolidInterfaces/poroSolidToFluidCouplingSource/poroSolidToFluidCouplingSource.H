/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::poroSolidToFluidCouplingSource

Group
    grpFvOptionsSources

Description
    Applies the kinematic source term from solid deformation onto the fluid
    mass conservation equation.

Author
    Denis Maier, BAW. All rights reserved.
    
See also
  - Foam::fvOption

SourceFiles
    poroSolidToFluidCouplingSource.C

\*---------------------------------------------------------------------------*/

#ifndef poroSolidToFluidCouplingSource_H
#define poroSolidToFluidCouplingSource_H

#include "Tuple2.H"
#include "cellSetOption.H"
#include "Enum.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declarations
class poroSolidToFluidCouplingSource;

Ostream& operator<<
(
    Ostream&,
    const poroSolidToFluidCouplingSource&
);


/*---------------------------------------------------------------------------*\
                     Class poroSolidToFluidCouplingSource Declaration
\*---------------------------------------------------------------------------*/

class poroSolidToFluidCouplingSource
:
    public fv::cellSetOption
{
public:
    // Public Enumerations


protected:

    // Protected Data
    word interfaceName_;
    // Protected Functions

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("poroSolidToFluidCouplingSource");


    // Constructors

        //- Construct from components
        poroSolidToFluidCouplingSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Access
        // Edit

        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            );

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
