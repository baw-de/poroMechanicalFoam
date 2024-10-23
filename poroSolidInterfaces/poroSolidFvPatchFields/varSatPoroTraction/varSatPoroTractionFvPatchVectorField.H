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
    varSatPoroTractionFvPatchVectorField

Description
    Fixed total traction boundary condition for partially saturated porous solvers.
    total traction is the sum of effective traction and (excess) pore water pressure.
    The contribution of the pore water pressure is controlled by a bishop function.

Author
    Denis Maier, BAW. All rights reserved.
    Héctor Mentenegro, BAW. All rights reserved.

    Based on the code by:
    Tian Tang, DTU.  
    Philip Cardiff, UCD. 

SourceFiles
    varSatPoroTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroTractionFvPatchVectorField_H
#define varSatPoroTractionFvPatchVectorField_H

#ifdef FOAMEXTEND
#include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "PatchFunction1.H"
#include "effectiveStressModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                 Class varSatPoroTractionFvPatchVectorField Declaration
    \*---------------------------------------------------------------------------*/

    class varSatPoroTractionFvPatchVectorField
        : public fixedGradientFvPatchVectorField
    {

        // Private Data

        //- total (true) or effective (false) traction/pressure
        const bool totalTraction_;

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Traction time series
        autoPtr<PatchFunction1<vector>> tractionSeries_;

        //- Pressure time series
        autoPtr<PatchFunction1<scalar>> pressureSeries_;

        //- effective stress model/definition
        autoPtr<effectiveStressModel> effectiveStressModelPtr_;
        //effectiveStressModel effectiveStressModelPtr_;

        //- name of chosen effectiveStressModel
        word effectiveStressModelName_;

        //- Second order accuracy correction
        const Switch secondOrder_;

        //- Limiter coefficient for snGrad correction
        const scalar limitCoeff_;

        //- Optional: under-relaxation factor
        const scalar relaxFac_;

        //- Optional: buoayancy is treated implicitly
        bool buoyancyImplicit_;

        effectiveStressModel& effectiveStressModelRef();

        void makeEffectiveStressModel();

    public:
        //- Runtime type information
        TypeName("varSatPoroTraction");

        // Constructors

        //- Construct from patch and internal field
        varSatPoroTractionFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &);

        //- Construct from patch, internal field and dictionary
        varSatPoroTractionFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  varSatPoroTractionFvPatchVectorField onto a new patch
        varSatPoroTractionFvPatchVectorField(
            const varSatPoroTractionFvPatchVectorField &,
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        varSatPoroTractionFvPatchVectorField(
            const varSatPoroTractionFvPatchVectorField &);

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>(
                new varSatPoroTractionFvPatchVectorField(*this));
        }

        //- Construct as copy setting internal field reference
        varSatPoroTractionFvPatchVectorField(
            const varSatPoroTractionFvPatchVectorField &,
            const DimensionedField<vector, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(
            const DimensionedField<vector, volMesh> &iF) const
        {
            return tmp<fvPatchVectorField>(
                new varSatPoroTractionFvPatchVectorField(*this, iF));
        }

        // Member functions

        // Access

        virtual const vectorField &traction() const
        {
            return traction_;
        }

        virtual vectorField &traction()
        {
            return traction_;
        }

        virtual const scalarField &pressure() const
        {
            return pressure_;
        }

        virtual scalarField &pressure()
        {
            return pressure_;
        }

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
        const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream &) const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
