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
    Foam::saturationLaws

Description
    A namespace for various saturationLaw implementations.

Class
    Foam::saturationLaw

Description
    An abstract base class for saturationLaws

    Constitutive relations for soils: S(h) and K(h)

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    saturationLaw.C
    newsaturationLaw.C

\*---------------------------------------------------------------------------*/

#ifndef saturationLaw_H
#define saturationLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                Class saturationLaw Declaration
\*---------------------------------------------------------------------------*/

    class saturationLaw
    {

    protected:
        // Protected data

        word name_;

        const volScalarField &pField_;

        // Private Member Functions

        //void makeDkbydp();

        //- Disallow copy construct
        saturationLaw(const saturationLaw &);

        //- Disallow default bitwise assignment
        void operator=(const saturationLaw &);

    public:
        //- Runtime type information
        TypeName("saturationLaw");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            saturationLaw,
            dictionary,
            (
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField),
            (name, poroHydraulicProperties, pField));

        // Selectors

        //- Return a reference to the selected soil model
        static autoPtr<saturationLaw> New(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField);

        // Constructors

        //- Construct from components
        saturationLaw(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField);

        // Destructor

        virtual ~saturationLaw()
        {
        }

        // Member Functions

        //- Read transportProperties dictionary

        virtual word name()
        {
            return name_;
        }

        // Access

        const fvMesh &mesh() const
        {
            return pField_.mesh();
        }

        const objectRegistry &db() const
        {
            return pField_.db();
        }

        virtual scalar pStar(const label cellI) const = 0;

        virtual scalar C(const scalar p, const label cellI) = 0;

        virtual scalar S(const scalar p, const label cellI) = 0;

        virtual scalar S
        (
            const scalar p,
            const label patchI,
            const label faceI
        ) = 0;

        virtual scalar kr(const scalar p, const label cellI) = 0;

        virtual scalar kr
        (
            const scalar p,
            const label patchI,
            const label faceI
        ) = 0;


        virtual void write() const = 0;

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
