/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::limitedHeadInfiltrationFvPatchScalarField

Description
    Foam::limitedHeadInfiltrationFvPatchScalarField

SourceFiles
    limitedHeadInfiltrationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef limitedHeadInfiltrationFvPatchScalarField_H
#define limitedHeadInfiltrationFvPatchScalarField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "mixedFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class limitedHeadInfiltrationFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class limitedHeadInfiltrationFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the z field
        word pName_;

        //- Name of the flux transporting the field
        word kEffname_;

        //- Name of the flux transporting the field
        word gradZname_;

        //- when was the last update?
        scalar updatedTime_;

        //- Total temperature
        scalarField pMax_;

        //- Pressure time series
        autoPtr<PatchFunction1<scalar>> fluxSeries_;

        //- is field pressure head or pressure
        Switch isHead_;

        //- where to look for fields when hydro-mechanically coupled
        word HMCoupled;
        
        //- recalculate weights implicitly (every iteration) or explicitly (every timestep)
        Switch	explicit_;
        

public:

    //- Runtime type information
    TypeName("limitedHeadInfiltration");


    // Constructors

        //- Construct from patch and internal field
        limitedHeadInfiltrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        limitedHeadInfiltrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  limitedHeadInfiltrationFvPatchScalarField onto a new patch
        limitedHeadInfiltrationFvPatchScalarField
        (
            const limitedHeadInfiltrationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        limitedHeadInfiltrationFvPatchScalarField
        (
            const limitedHeadInfiltrationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new limitedHeadInfiltrationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        limitedHeadInfiltrationFvPatchScalarField
        (
            const limitedHeadInfiltrationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new limitedHeadInfiltrationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the total pressure
            const scalarField& pMax() const
            {
                return pMax_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& pMax()
            {
                return pMax_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
