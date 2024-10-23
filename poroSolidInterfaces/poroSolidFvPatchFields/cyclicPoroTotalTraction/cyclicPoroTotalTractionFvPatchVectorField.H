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
    cyclicPoroTotalTractionFvPatchVectorField

Description
    Fixed traction boundary condition for solid solvers.

Author
    Denis Maier, BAW.
    Based on the work by Tian Tang, DTU.

SourceFiles
    cyclicPoroTotalTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPoroTotalTractionFvPatchVectorField_H
#define cyclicPoroTotalTractionFvPatchVectorField_H

#ifdef FOAMEXTEND
#include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
             Class cyclicPoroTotalTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

    class cyclicPoroTotalTractionFvPatchVectorField
        : public fixedGradientFvPatchVectorField
    {

        // Private Data

        //- Traction
        vectorField tractionAmp_; // traction amplitude

        //- Pressure
        scalarField pressureAmp_; // pressure amplitude

        //- Traction
        vectorField tractionBase_; // Base traction

        //- Pressure
        scalarField pressureBase_; // Base pressure

        //- Period
        scalar T_;

        //- Frequency
        scalar f_;

        //- Frequency
        scalar shift_;

        //- Second order accuracy correction
        const Switch secondOrder_;

        //- Limiter coefficient for snGrad correction
        const scalar limitCoeff_;

        //- Optional: under-relaxation factor
        const scalar relaxFac_;

    public:
        //- Runtime type information
        TypeName("cyclicPoroTotalTraction");

        // Constructors

        //- Construct from patch and internal field
        cyclicPoroTotalTractionFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &);

        //- Construct from patch, internal field and dictionary
        cyclicPoroTotalTractionFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  cyclicPoroTotalTractionFvPatchVectorField onto a new patch
        cyclicPoroTotalTractionFvPatchVectorField(
            const cyclicPoroTotalTractionFvPatchVectorField &,
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        cyclicPoroTotalTractionFvPatchVectorField(
            const cyclicPoroTotalTractionFvPatchVectorField &);

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>(
                new cyclicPoroTotalTractionFvPatchVectorField(*this));
        }

        //- Construct as copy setting internal field reference
        cyclicPoroTotalTractionFvPatchVectorField(
            const cyclicPoroTotalTractionFvPatchVectorField &,
            const DimensionedField<vector, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(
            const DimensionedField<vector, volMesh> &iF) const
        {
            return tmp<fvPatchVectorField>(
                new cyclicPoroTotalTractionFvPatchVectorField(*this, iF));
        }

        // Member functions

        // Access

        virtual const vectorField &tractionAmp() const
        {
            return tractionAmp_;
        }

        virtual vectorField &tractionAmp()
        {
            return tractionAmp_;
        }

        virtual const scalarField &pressureAmp() const
        {
            return pressureAmp_;
        }

        virtual scalarField &pressureAmp()
        {
            return pressureAmp_;
        }

        virtual const vectorField &tractionBase() const
        {
            return tractionBase_;
        }

        virtual vectorField &tractionBase()
        {
            return tractionBase_;
        }

        virtual const scalarField &pressureBase() const
        {
            return pressureBase_;
        }

        virtual scalarField &pressureBase()
        {
            return pressureBase_;
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
