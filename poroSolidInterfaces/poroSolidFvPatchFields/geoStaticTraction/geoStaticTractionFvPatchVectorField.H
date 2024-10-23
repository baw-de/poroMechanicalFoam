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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    geoStaticTractionFvPatchVectorField

Description
    Fixed traction boundary condition for solid solvers where the traction
    is generated from the solid density.

Author
    Denis Maier, BAW.
SourceFiles
    geoStaticTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef geoStaticTractionFvPatchVectorField_H
#define geoStaticTractionFvPatchVectorField_H

#ifdef FOAMEXTEND
#include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "poroTractionFvPatchVectorField.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
             Class geoStaticTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

    class geoStaticTractionFvPatchVectorField
        : public poroTractionFvPatchVectorField
    {

        // Private Data

        //- z GOK
        scalar K0_;

        //- spezific weight
        scalar gamma_;

        //- Location of ground level
        scalar GOK_;

        //- excavation elevation time series
        interpolationTable<scalar> zSeries_;

    public:
        //- Runtime type information
        TypeName("geoStaticTraction");

        // Constructors

        //- Construct from patch and internal field
        geoStaticTractionFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &);

        //- Construct from patch, internal field and dictionary
        geoStaticTractionFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  geoStaticTractionFvPatchVectorField onto a new patch
        geoStaticTractionFvPatchVectorField(
            const geoStaticTractionFvPatchVectorField &,
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        geoStaticTractionFvPatchVectorField(
            const geoStaticTractionFvPatchVectorField &);

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>(
                new geoStaticTractionFvPatchVectorField(*this));
        }

        //- Construct as copy setting internal field reference
        geoStaticTractionFvPatchVectorField(
            const geoStaticTractionFvPatchVectorField &,
            const DimensionedField<vector, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(
            const DimensionedField<vector, volMesh> &iF) const
        {
            return tmp<fvPatchVectorField>(
                new geoStaticTractionFvPatchVectorField(*this, iF));
        }

        // Member functions

        // Access

        // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(
            const fvPatchFieldMapper &);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap(
            const fvPatchVectorField &,
            const labelList &);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate(
#ifdef OPENFOAMESIORFOUNDATION
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );

        //- Write
        virtual void write(Ostream &) const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
