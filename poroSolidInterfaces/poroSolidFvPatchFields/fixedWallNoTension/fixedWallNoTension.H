/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    fixedWallNoTension

Description
    Same as fixedDisplacement, but changes to zero traction when stress would be 
    induced by the fixedDisplacement boundary condition.

SourceFiles
    fixedWallNoTension.C

Author
    Denis Maier, BAW. All rights reserved.
    Based on the code by Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef fixedWallNoTension_H
#define fixedWallNoTension_H

#include "foamTime.H"
#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

    class fixedWallNoTension
        : public solidDirectionMixedFvPatchVectorField
    {
        // Private data

        //- Constant displacement field

    public:
        //- Runtime type information
        TypeName("fixedWallNoTension");

        // Constructors

        //- Construct from patch and internal field
        fixedWallNoTension(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &);

        //- Construct from patch, internal field and dictionary
        fixedWallNoTension(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        // fixedWallNoTension
        //  onto a new patch
        fixedWallNoTension(
            const fixedWallNoTension &,
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>(
                new fixedWallNoTension(*this));
        }

        //- Construct as copy setting internal field reference
        fixedWallNoTension(
            const fixedWallNoTension &,
            const DimensionedField<vector, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone(
            const DimensionedField<vector, volMesh> &iF) const
        {
            return tmp<fvPatchField<vector>>(
                new fixedWallNoTension(*this, iF));
        }

        // Destructor

        virtual ~fixedWallNoTension()
        {
        }

        // Member functions

        // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(
            const fvPatchFieldMapper &);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap(
            const fvPatchField<vector> &,
            const labelList &);

        // Access functions

        // Evaluation functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream &) const;

        // Member operators
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
