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

\*---------------------------------------------------------------------------*/

#include "geoStaticTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    geoStaticTractionFvPatchVectorField::
        geoStaticTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF)
        : poroTractionFvPatchVectorField(p, iF),
          K0_(0.0),
          gamma_(0.0),
          GOK_(0.0),
          zSeries_()
    {
        fvPatchVectorField::operator=(patchInternalField());
        gradient() = vector::zero;
    }

    geoStaticTractionFvPatchVectorField::
        geoStaticTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const dictionary &dict)
        : poroTractionFvPatchVectorField(p, iF),
          K0_(readScalar(dict.lookup("K0"))),
          gamma_(readScalar(dict.lookup("gamma"))),
          GOK_(0.0),
          zSeries_()
    {
        Info << "Creating " << type() << " boundary condition" << endl;
        if (dict.found("GroundLevelSeries"))
        {
            zSeries_ = interpolationTable<scalar>(dict.subDict("GroundLevelSeries"));
        }
        else
        {
            GOK_ = readScalar(dict.lookup("GroundLevel"));
        }

        if (dict.found("gradient"))
        {
            gradient() = vectorField("gradient", dict, p.size());
        }
        else
        {
            gradient() = vector::zero;
        }

        if (dict.found("value"))
        {
            Field<vector>::operator=(vectorField("value", dict, p.size()));
        }
        else
        {
            fvPatchVectorField::operator=(patchInternalField());
        }
    }

    geoStaticTractionFvPatchVectorField::
        geoStaticTractionFvPatchVectorField(
            const geoStaticTractionFvPatchVectorField &stpvf,
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const fvPatchFieldMapper &mapper)
        : poroTractionFvPatchVectorField(stpvf, p, iF, mapper),
          K0_(stpvf.K0_),
          gamma_(stpvf.gamma_),
          GOK_(stpvf.GOK_),
          zSeries_(stpvf.zSeries_)
    {
    }

    geoStaticTractionFvPatchVectorField::
        geoStaticTractionFvPatchVectorField(
            const geoStaticTractionFvPatchVectorField &stpvf)
        : poroTractionFvPatchVectorField(stpvf),
          K0_(stpvf.K0_),
          gamma_(stpvf.gamma_),
          GOK_(stpvf.GOK_),
          zSeries_(stpvf.zSeries_)
    {
    }

    geoStaticTractionFvPatchVectorField::
        geoStaticTractionFvPatchVectorField(
            const geoStaticTractionFvPatchVectorField &stpvf,
            const DimensionedField<vector, volMesh> &iF)
        : poroTractionFvPatchVectorField(stpvf, iF),
          K0_(stpvf.K0_),
          gamma_(stpvf.gamma_),
          GOK_(stpvf.GOK_),
          zSeries_(stpvf.zSeries_)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void geoStaticTractionFvPatchVectorField::autoMap(
        const fvPatchFieldMapper &m)
    {
        poroTractionFvPatchVectorField::autoMap(m);
    }

    // Reverse-map the given fvPatchField onto this fvPatchField
    void geoStaticTractionFvPatchVectorField::rmap(
        const fvPatchVectorField &ptf,
        const labelList &addr)
    {
        poroTractionFvPatchVectorField::rmap(ptf, addr);

        const geoStaticTractionFvPatchVectorField &dmptf =
            refCast<const geoStaticTractionFvPatchVectorField>(ptf);
    }

    // Update the coefficients associated with the patch field
    void geoStaticTractionFvPatchVectorField::updateCoeffs()
    {
        scalarField z(patch().size(), 0.0);
        if (updated())
        {
            return;
        }

        /*if (solMod.name() == "sixDoFRigidSolid")
        {
            const fvsPatchField<vector> &C =
                lookupPatchField<surfaceVectorField, vector>("faceLocation");
            z = C.component(vector::Z);
        }
        else
        {*/
            const vectorField &C = patch().Cf();
            z = C.component(vector::Z);
        //}

        if (zSeries_.size())
        {
            GOK_ = zSeries_(db().time().timeOutputValue());
        }

        // when mesh is moving this is gonna be implemented instead: scalarField z(patch().Cf().component(vector::Z));
        vectorField traction(patch().size(), vector::zero);
        scalarField pressure(patch().size(), 0.0);

        // symmTensorField sigmaStatic(patch().size(),symmTensor::zero);
        vectorField n = patch().nf();

        forAll(patch(), iFace)
        {
            if (z[iFace] < GOK_)
            {
                // sigmaStatic[iFace].zz() = -mag(gamma_) * (GOK_ - z[iFace]);
                // sigmaStatic[iFace].xx() = -K0_ * mag(gamma_) * (GOK_ - z[iFace]);
                // sigmaStatic[iFace].yy() = -K0_ * mag(gamma_) * (GOK_ - z[iFace]);
                pressure[iFace] = -K0_ * mag(gamma_) * (GOK_ - z[iFace]);
            }
            else
            {
                // sigmaStatic[iFace].zz() = 0.0;
                // sigmaStatic[iFace].xx() = 0.0;
                // sigmaStatic[iFace].yy() = 0.0;
                pressure[iFace] = 0.0;
            }
        }
        // Lookup the solidModel object
        // traction =  sigmaStatic & n;

        // Set surface-normal gradient on the patch corresponding to the desired
        // traction
        poroTractionFvPatchVectorField::pressure() = pressure;
        poroTractionFvPatchVectorField::updateCoeffs();
    }

    void geoStaticTractionFvPatchVectorField::evaluate(
        const Pstream::commsTypes commsType)
    {
        poroTractionFvPatchVectorField::evaluate();
    }

    void geoStaticTractionFvPatchVectorField::write(Ostream &os) const
    {

        if (zSeries_.size())
        {
            os.writeKeyword("GOKSeries") << nl;
            os << token::BEGIN_BLOCK << nl;
            zSeries_.write(os);
            os << token::END_BLOCK << nl;
        }
        else
        {
            os << tab << "GroundLevel" << tab << GOK_ << ";" << endl;
        }

        os.writeKeyword("K0")
            << K0_ << token::END_STATEMENT << nl;
        os.writeKeyword("gamma")
            << gamma_ << token::END_STATEMENT << nl;
        poroTractionFvPatchVectorField::write(os);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField(fvPatchVectorField, geoStaticTractionFvPatchVectorField);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
