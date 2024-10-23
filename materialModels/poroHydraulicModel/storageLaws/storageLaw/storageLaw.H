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

Namespace
    Foam::storageLaws

Description
    A namespace for various storageLaw implementations.

Class
    Foam::storageLaw

Description
    An abstract base class for storageLaws

    Constitutive relations for soils: Ss(h)
    
Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    storageLaw.C
    newstorageLaw.C

\*---------------------------------------------------------------------------*/

#ifndef storageLaw_H
#define storageLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                Class storageLaw Declaration
\*---------------------------------------------------------------------------*/

    class storageLaw
    {
        // Protected data
        // name of the storage law
        word name_;
        dictionary &poroHydraulicProperties_;

    protected:

        //- Pressure field used by solver
        const volScalarField &pField_;

        //- solid compressibility
        autoPtr<volScalarField> Cs_;

        bool writeSs_;

        // Private Member Functions
        void makeCs();

        //- Disallow copy construct
        storageLaw(const storageLaw &);

        //- Disallow default bitwise assignment
        void operator=(const storageLaw &);

    public:
        //- Runtime type information
        TypeName("storageLaw");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            storageLaw,
            dictionary,
            (
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField),
            (name, poroHydraulicProperties, pField));

        // Selectors

        //- Return a reference to the selected soil model
        static autoPtr<storageLaw> New(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField);

        // Constructors

        //- Construct from components
        storageLaw(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField);

        // Destructor

        virtual ~storageLaw()
        {
        }

        // Member Functions

        const objectRegistry &db() const
        {
            return pField_.db();
        }

        const fvMesh &mesh() const
        {
            return pField_.mesh();
        }

        //- Read transportProperties dictionary

        virtual word name()
        {
            return name_;
        }

        virtual scalar Ss
        (
            const scalar n,
            const scalar p,
            const label cellI
        ) = 0;

        virtual bool updatesSs() const 
        {
            return true;
        }

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
