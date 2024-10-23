/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Namespace
    Foam::incompressible::scalarDiskReaders

Description
    Namespace for incompressible turbulence turbulence models.

Class
    Foam::incompressible::scalarDiskReader

Description
    Abstract base class for incompressible turbulence models
    (RAS, LES and laminar).

Authors
    Heavily inspired by code of Philip Cardiff, UCD. All rights reserved.
    Denis Maier, BAW, packaging into seperate class format.

SourceFiles
    scalarDiskReader.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef scalarDiskReader_H
#define scalarDiskReader_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "meshToMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                      Class scalarDiskReader Declaration
\*---------------------------------------------------------------------------*/

    class scalarDiskReader
    {

        // Private data

        //- Name of the field to be read from file
        const word scalarName_;

        //- Reference to mesh
        const fvMesh &mesh_;

        //- Reference to baseMesh
        const fvMesh &baseMesh_;

        //- meshToMesh mapping
        autoPtr<meshToMesh> scalarToSolidMapping_;

        //- dictionary which might contain meshToMesh mapping data
        dictionary dict_;

        //- field read from file (on solidSubMesh)
        autoPtr<volScalarField> FieldPtr_;

        //- field read from file (on baseMesh)
        autoPtr<volScalarField> fieldBaseMeshPtr_;

        //- flag indicating if scalar value was read from file
        bool fieldWasReadFromDisk_;

        //- Path to case directory where scalar data is held
        const word fieldCaseDir_;

        //- new runtime to read data from
        Time fieldRunTime;

        //- mesh on which the original data is held
        fvMesh fieldMesh_;

        //- time index to indicate when the last update occured
        label curTimeIndex_;

        // Private Member Functions
        void makeMeshMapping();

        // return meshToMesh mapping
        meshToMesh &scalarToSolidMapping();

        //- Disallow default bitwise copy construct
        scalarDiskReader(const scalarDiskReader &);

        //- Disallow default bitwise assignment
        void operator=(const scalarDiskReader &);

    public:

        // Constructors

        //- Construct from components
        scalarDiskReader(
        const word scalarName,
        const fvMesh &mesh,
        const fvMesh &baseMesh,
        const dictionary dict);

        // Destructor

        virtual ~scalarDiskReader()
        {
          scalarToSolidMapping_.clear();
          FieldPtr_.clear();
          fieldBaseMeshPtr_.clear();
        }

        // Member Functions
        bool readField();

        //- Return the read field
        tmp<volScalarField> field()
        {
            // Only attempt to read the T field from disk once per time-step
            if (curTimeIndex_ != mesh_.time().timeIndex())
            {
                curTimeIndex_ = mesh_.time().timeIndex();
                fieldWasReadFromDisk_ = readField();
            }
            return FieldPtr_();
        }

        //- Return the read field
        tmp<volScalarField> field() const
        {
            if(!FieldPtr_.valid())
            {
                FatalErrorInFunction() << "ScalarField " << scalarName_
                                       << " has not yet been read from disk!"
                                       << endl;
            }
            return FieldPtr_();
        }

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
