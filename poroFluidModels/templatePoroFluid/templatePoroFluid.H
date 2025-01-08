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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::poroFluidModels::[template]

Description
    Template for implementing a groundwater flow solver using pressure deviation 
    from initial hydrostatic conditions (p_rgh).
    
    This solver works with pressure units rather than water column height, making
    it equivalent to the groundwater potential form but in different units.

    Key Features:
    - Runtime-selectable storage and hydraulic conductivity functions via poroHydraulicModel
    - Handles both saturated and unsaturated flow conditions
    - Integrated with OpenFOAM's infrastructure for parallel computing and mesh handling

    Implementation Notes:
    1. Replace '[template]' with your desired model name throughout this file
    2. The solver uses p_rgh (dynamic pressure) as the primary variable

Required Files:
    - [template].C : Implementation file

Author
    Denis Maier, BAW. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef [template]_H
#define [template]_H

#include "poroFluidModel.H"
#include "poroHydraulicModel.H"

namespace Foam
{
namespace poroFluidModels
{

/*---------------------------------------------------------------------------*\
                        Class [template] Declaration
\*---------------------------------------------------------------------------*/

class [template]
:
    public poroFluidModel
{
    // Private Data Members
    
    //- Pointer to the hydraulic transport model
    //  Handles conductivity and storage properties
    mutable autoPtr<poroHydraulicModel> poroHydPtr_;

    //- Storage coefficient field [1/m]
    //  Represents the pressure reaction to in/outflow and compression per unit volume 
    volScalarField Ss_;

    //- Surface field of effective hydraulic conductivity divided by water specific weight
    //  Pre-computed for efficiency in flux calculations
    surfaceScalarField kbyGammaf_;

    // Private Member Functions

    //- Access to hydraulic model
    poroHydraulicModel& poroHydraulic();

    //- Disallow default bitwise copy construct
    [template](const [template]&);

    //- Disallow default bitwise assignment
    void operator=(const [template]&);

public:
    //- Runtime type information for dynamic typing
    TypeName("[template]");

    // Constructors

    //- Construct from components
    //  @param runTime Reference to OpenFOAM time object
    //  @param region Mesh region name
    //  @param sharedMesh Flag indicating if mesh is shared
    [template]
    (
        Time& runTime,
        const word& region = dynamicFvMesh::defaultRegion,
        const bool sharedMesh = false
    );

    //- Destructor
    virtual ~[template]()
    {}

    // Member Functions

    // Access

    //- Access to pressure field (convenience wrapper)
    virtual volScalarField& pField()
    {
        return p_rgh();
    }
    
    //- Const access to pressure field
    virtual const volScalarField& pField() const
    {
        return p_rgh();
    }

    //- Return initial porosity field
    virtual tmp<volScalarField> n0()
    {
        return poroHydraulic().n0();
    }

    //- Const access to hydraulic model
    const poroHydraulicModel& poroHydraulic() const
    {
        return poroHydPtr_();
    }

    // Evolution and Control

    //- Solve for one time step
    //  @return true if solution converged
    virtual bool evolve();

    //- Calculate new time step size based on solution behavior
    //  @return Suggested next time step size
    virtual scalar newDeltaT();

    //- Return material residual (always 0.0 for linear equations)
    virtual scalar materialResidual()
    {
        return 0.0;
    }

    //- Return number of outer iterations
    virtual label outerIters()
    {
        return iterCtrl().index();
    }

    //- Write fields to disk
    virtual void writeFields(const Time& runTime);

    //- Clean up at end of simulation
    virtual void end();
};

} // End namespace poroFluidModels
} // End namespace Foam

#endif

// ************************************************************************* //