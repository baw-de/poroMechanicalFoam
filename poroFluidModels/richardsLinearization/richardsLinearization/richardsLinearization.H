/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Namespace
    Foam::richardsLinearizations

Description
    A namespace for various richardsLinearization implementations.

Class
    Foam::richardsLinearization

Description
    An abstract base class for linearization schemes for 
    Richard's equation.

Author
    Denis Maier, BAW.   All rights reserved.
    
SourceFiles
    richardsLinearization.C
    newrichardsLinearization.C

\*---------------------------------------------------------------------------*/

#ifndef richardsLinearization_H
#define richardsLinearization_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "volFields.H"
#include "varSatPoroHydraulicModel.H"
#ifdef OPENFOAMESIORFOUNDATION
#include "SolverPerformance.H"
#endif
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                Class richardsLinearization Declaration
\*---------------------------------------------------------------------------*/

    class richardsLinearization
    {

    protected:
        // Protected data

        word name_;

        varSatPoroHydraulicModel &poroHydraulic_;

        volScalarField& S_;

        volScalarField C_;

        const fvMesh &mesh_;

        // Private Member Functions

        //- Disallow copy construct
        richardsLinearization(const richardsLinearization &);

        //- Disallow default bitwise assignment
        void operator=(const richardsLinearization &);

        volScalarField& S()
        {
            return S_;
        }

    public:
        //- runTime type information
        TypeName("richardsLinearization");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            richardsLinearization,
            dictionary,
            (
                const word &name,
                varSatPoroHydraulicModel &poroHydraulic,
                dictionary &poroFluidProperties,
                volScalarField& S),
            (name, poroHydraulic, poroFluidProperties, S));

        // Selectors

        //- Return a reference to the selected poroFluid model
        static autoPtr<richardsLinearization> New(
            const word &name,
            varSatPoroHydraulicModel &poroHydraulic,
            dictionary &poroFluidProperties,
            volScalarField& S);

        // Constructors

        //- Construct from components
        richardsLinearization(
            const word &name,
            varSatPoroHydraulicModel &poroHydraulic,
            dictionary &poroFluidProperties,
            volScalarField& S);

        // Destructor

        virtual ~richardsLinearization()
        {
        }

        // Member Functions

        virtual const word name() const {
            return name_;
        }

        const fvMesh &mesh() const
        {
            return mesh_;
        }

        const volScalarField& S() const
        {
            return S_;
        }

        volScalarField& C()
        {
            return C_;
        }

        const volScalarField& C() const
        {
            return C_;
        }

        varSatPoroHydraulicModel &poroHydraulic()
        {
            return poroHydraulic_;
        }

        const varSatPoroHydraulicModel &poroHydraulic() const
        {
            return poroHydraulic_;
        }

        virtual void initalize(volScalarField &totalP,volScalarField &pField) = 0;

        virtual bool checkConvergedAndUpdate(volScalarField &totalP, volScalarField &pField) = 0;

        virtual tmp<fvScalarMatrix> ddtS(const volScalarField &S, volScalarField &pField) = 0;
        
        //tmp<fvScalarMatrix> ddpk(const surfaceScalarField &kField, volScalarField &pField);

        scalar L2Norm(const scalarField &x) const;

        /*scalar MassBalance(const volScalarField &pField, const volScalarField &n, const volScalarField &S, const volScalarField &flux) const;

        tmp<scalarField> volMassChange(const volScalarField &theta) const;

        tmp<scalarField> volMiscStorage(const volScalarField &pField, const volScalarField &n) const;

        tmp<volScalarField> massFlux(const surfaceScalarField &k_eff, const volScalarField &pField) const;*/
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
