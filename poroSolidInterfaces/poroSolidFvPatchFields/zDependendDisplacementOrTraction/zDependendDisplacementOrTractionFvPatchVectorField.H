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
    zDependendDisplacementOrTractionFvPatchVectorField

Description
    Time varying Displacment an traction boundary condition, that devides the 
    patch into two zones. A fixed displacment and a traction zone. They are 
    controlled by a function z(t), where z is compared to the z-coordinate of the
    face.

SourceFiles
    zDependendDisplacementOrTractionFvPatchVectorField.C

Author
    Denis Maier, BAW.All rights reserved.
    Based on the code by:
    Philip Cardiff, UCD. 
\*---------------------------------------------------------------------------*/

#ifndef zDependendDisplacementOrTractionFvPatchVectorField_H
#define zDependendDisplacementOrTractionFvPatchVectorField_H

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

class zDependendDisplacementOrTractionFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Constant displacement field
        vectorField constantDisplacement_;

        //- Constant traction field
        vectorField constantTraction_;

        //- Constant pressure field
        scalarField constantPressure_;

        //- Constant boundary elevation
        scalar z0_;

        //- Switch that allows for specifiing which direction the split has
        Switch splitDir_;

        //- Optional: displacement time series
        interpolationTable<vector> displacementSeries_;

        //- Optional: traction time series
        interpolationTable<vector> tractionSeries_;

        //- Optional: pressure time series
        interpolationTable<scalar> pressureSeries_;

        //- Optional: boundary elevation time series
        interpolationTable<scalar> zSeries_;

public:

    //- Runtime type information
    TypeName("zDependendDisplacementOrTraction");


    // Constructors

        //- Construct from patch and internal field
        zDependendDisplacementOrTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zDependendDisplacementOrTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // zDependendDisplacementOrTractionFvPatchVectorField
        //  onto a new patch
        zDependendDisplacementOrTractionFvPatchVectorField
        (
            const zDependendDisplacementOrTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new zDependendDisplacementOrTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zDependendDisplacementOrTractionFvPatchVectorField
        (
            const zDependendDisplacementOrTractionFvPatchVectorField&,
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
                new zDependendDisplacementOrTractionFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~zDependendDisplacementOrTractionFvPatchVectorField()
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

            // Const access to displacement series
            const interpolationTable<vector>& displacementSeries() const
            {
                return displacementSeries_;
            }

            // Const access to traction series
            const interpolationTable<vector>& tractionSeries() const
            {
                return tractionSeries_;
            }

            // Const access to pressure series
            const interpolationTable<scalar>& pressureSeries() const
            {
                return pressureSeries_;
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
