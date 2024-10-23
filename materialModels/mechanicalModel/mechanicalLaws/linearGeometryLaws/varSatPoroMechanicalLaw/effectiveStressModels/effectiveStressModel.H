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
    Foam::effectiveStressModels

Description
    A namespace for various effectiveStressModel implementations.

Class
    Foam::effectiveStressModel

Description
    An abstract base class for effectiveStressModels

    Constitutive relations for soils: Ss(h)

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    effectiveStressModel.C
    neweffectiveStressModel.C

\*---------------------------------------------------------------------------*/

#ifndef effectiveStressModel_H
#define effectiveStressModel_H

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
                Class effectiveStressModel Declaration
\*---------------------------------------------------------------------------*/

    class effectiveStressModel
    {

        // Protected data
        // name of the storage law
        word name_;

        const fvMesh &mesh_;

    protected:

        autoPtr<volScalarField> chi_;

        void makeChi();

        const fvMesh &mesh() const
        {
            return mesh_;
        }

        const dictionary &effectiveStressModelDict_;

        //- Disallow copy construct
        effectiveStressModel(const effectiveStressModel &);

        //- Disallow default bitwise assignment
        void operator=(const effectiveStressModel &);

    public:
        //- Runtime type information
        TypeName("effectiveStressModel");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            effectiveStressModel,
            dictionary,
            (
                const dictionary &effectiveStressModelDict,
                const word &name,
                const fvMesh &mesh),
            (effectiveStressModelDict, name, mesh));

        // Selectors

        //- Return a reference to the selected soil model
        static autoPtr<effectiveStressModel> New(
            const dictionary &effectiveStressModelDict,
            const word &name,
            const fvMesh &mesh);

        // Constructors

        //- Construct from components
        effectiveStressModel(
            const dictionary &effectiveStressModelDict,
            const word &name,
            const fvMesh &mesh);

        // Destructor

        virtual ~effectiveStressModel()
        {
        }

        // Member Functions

        //- Read transportProperties dictionary

        virtual word name()
        {
            return name_;
        }

        virtual tmp<volScalarField> chi()
        {
            if(!chi_.valid())
            {
                makeChi();
            }
            return chi_();
        }

        virtual tmp<volScalarField> chi(const volScalarField &n, const volScalarField &S, const volScalarField &p);

        virtual tmp<scalarField> chi(const scalarField &n, const scalarField &S, const scalarField &p) = 0;

        virtual void write() const {};
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
