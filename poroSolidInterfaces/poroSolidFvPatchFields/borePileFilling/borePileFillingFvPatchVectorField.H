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
    borePileFillingFvPatchVectorField

Description
    This boundary condition allows a total pressure or zero displacement 
    to be specified in different zones, depending on whether a function z(t) 
    lies under the z-coordinate of the face or not. 
    This simulates a bore pile that is filled with concrete and is being
    pulled out of the ground.

SourceFiles
    borePileFillingFvPatchVectorField.C

Author
    Denis Maier, BAW. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef borePileFillingFvPatchVectorField_H
#define borePileFillingFvPatchVectorField_H

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

class borePileFillingFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data
        //- Constant pile displacement = 0
        vectorField constantDisplacement_;

        //- Constant pile elevation
        scalar z0_;

        //- Constant concrete head
        scalar h_c_;

        //- gamma of concrete
        scalar gamma_c_;

        //- Switch that allows for specifiing which direction the split has
        Switch splitDir_;

        //- Optional: concrete head time series
        interpolationTable<scalar> hcSeries_;

        //- Optional: boundary elevation time series
        interpolationTable<scalar> zSeries_;

public:

    //- Runtime type information
    TypeName("borePileFilling");


    // Constructors

        //- Construct from patch and internal field
        borePileFillingFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        borePileFillingFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // borePileFillingFvPatchVectorField
        //  onto a new patch
        borePileFillingFvPatchVectorField
        (
            const borePileFillingFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new borePileFillingFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        borePileFillingFvPatchVectorField
        (
            const borePileFillingFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new borePileFillingFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~borePileFillingFvPatchVectorField()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );

        // Access functions


            // Const access to pressure series
            const interpolationTable<scalar>& hcSeries() const
            {
                return hcSeries_;
            }

            // Const access to z series
            const interpolationTable<scalar>& zSeries() const
            {
                return zSeries_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
