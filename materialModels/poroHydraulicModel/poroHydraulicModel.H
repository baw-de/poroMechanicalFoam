/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::poroHydraulicModel

Description
    The poroHydraulic model class takes care of reading the poroHydraulic properties
    and creating poroHydraulic laws for storage and hydraulic conductivity.

SourceFiles
    poroHydraulicModel.C

Author
    Denis Maier, BAW. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef poroHydraulicModel_H
#define poroHydraulicModel_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "saturationLaw.H"
#include "storageLaw.H"
#include "conductivityModel.H"
#include "Switch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
  /*---------------------------------------------------------------------------*\
                        Class poroHydraulicModel Declaration
\*---------------------------------------------------------------------------*/

  class poroHydraulicModel
      : public regIOobject  // so it can be registered, and found in time()
  {
    // Private data

    //- Dictrionary for porohydraulic coefficients
    IOdictionary poroHydraulicProperties_;

    //- Pressure field used by solver
    const volScalarField& pField_;

    //- law describing pore pressure storage term (e.g. from compressibility)
    mutable PtrList<storageLaw> storageLawPtr_;

    //- laws describing change of hydraulic conductivity changes (e.g. from porosity changes)
    mutable PtrList<conductivityModel> conductivityModelPtr_;

    //- Density of poreFluid
    const dimensionedScalar rho_;

    //- Specific weigth of poreFluid
    const UniformDimensionedField<vector> gamma_;

    //- reference level and derived gravitational potential fields
    //- reference water table (p_rgh is in reference to 
    //  the hydrostatic pressure generated by a watercolumn 
    //  of this height)
    const dimensionedScalar href_;
    
    //- Geodesic height  
    const volScalarField z_;

    //- hydrostatic pressure generated by a watercolumn of
    //  href_ height
    const volScalarField p_Hyd_;

    bool updatesSs_;
    bool updatesK_;

    //- Disallow default bitwise copy construct
    poroHydraulicModel(const poroHydraulicModel &);

    //- Disallow default bitwise assignment
    void operator=(const poroHydraulicModel &);

  protected:

    dictionary &poroHydraulicProperties()
    {
      return poroHydraulicProperties_;
    }

  public:
    //- Runtime type information
    TypeName("poroHydraulicModel");

    // Constructors
    //- Construct from an fvMesh
    poroHydraulicModel
    (
      const volScalarField& pField,
      const dimensionedVector& gravity
    );

    //- Destructor
    ~poroHydraulicModel();

    // Member Functions

    // Access

    const dictionary& poroHydraulicProperties() const
    {
      return poroHydraulicProperties_;
    }

    const volScalarField& pField() const
    {
      return pField_;
    }

    tmp<volScalarField> n0() const;

    //- Return const access to the mesh
    const fvMesh& mesh() const 
    {
        return pField_.mesh();
    }

    //- Return initial density
    const dimensionedScalar rho() const
    {
      return rho_;
    }

    const dimensionedScalar magGamma() const
    {
      return mag(gamma_);
    }

    const dimensionedVector gamma() const
    {
      return gamma_;
    }

    const dimensionedScalar href() const
    {
      return href_;
    }

    const volScalarField &z() const
    {
      return z_;
    }

    const volScalarField &p_Hyd() const
    {
      return p_Hyd_;
    }

    virtual tmp<volScalarField> Ss(const volScalarField& n, const volScalarField &p);

    const tmp<volScalarField> k() const;

    virtual const tmp<surfaceScalarField> kf() const;

    bool updatesSs()
    {
      return updatesSs_;
    }

    bool updatesK()
    {
      return updatesK_;
    }

    //- Write function must be defined for regIOobjects
    virtual bool writeData(Ostream&) const
    {
      return false;
    }

    //- Write function must be defined for regIOobjects
    virtual void write()
    {
      poroHydraulicProperties_.IOobject::rename(poroHydraulicProperties_.IOobject::name()+".withDefaultValues");
      poroHydraulicProperties_.regIOobject::write();
    }

    virtual void writeCoeffFields() const
    {
      // notImplemented;
    }
  };

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
