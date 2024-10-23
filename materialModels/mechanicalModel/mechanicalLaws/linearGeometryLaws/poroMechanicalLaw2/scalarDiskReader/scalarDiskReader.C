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

\*---------------------------------------------------------------------------*/

#include "scalarDiskReader.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * *  * * Private Member Functions  * * * * * * * * * * * * //

    void Foam::scalarDiskReader::makeMeshMapping()
    {

        const word mapMethodName(dict_.getOrDefault<word>("mapMethod","direct"));
        if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
        {
            FatalErrorInFunction
                << "unknown map method "
                << mapMethodName << nl
                << "Available methods include: "
                << meshToMesh::interpolationMethodNames_
                << exit(FatalError);
        }

        meshToMesh::interpolationMethod mapMethod
        (
            meshToMesh::interpolationMethodNames_[mapMethodName]
        );

        Switch consistent(
            dict_.getOrDefault<Switch>(
                "consistent", true));
            if(!consistent)
            {
                if(mapMethodName=="direct")
                {
                    WarningInFunction
                        << "direct mapping selected for inconsitent region meshes, changing to imMapNearest"
                        << endl;
                    mapMethod = meshToMesh::interpolationMethod::imMapNearest;
                }
            HashTable<word> patchMap;
            wordList cuttingPatches;

            dict_.readEntry("patchMap", patchMap);
            dict_.readEntry("cuttingPatches", cuttingPatches);
            scalarToSolidMapping_.reset
                    (
                        new meshToMesh
                        (
                            fieldMesh_,
                            baseMesh_,
                            mapMethod,
                            patchMap,
                            cuttingPatches
                        )
                    );
            }
            else
            {
            scalarToSolidMapping_.reset
                    (
                        new meshToMesh
                        (
                            fieldMesh_,
                            baseMesh_,
                            mapMethod
                        )
                    );
            }
    }

    Foam::meshToMesh &Foam::scalarDiskReader::scalarToSolidMapping()
    {
        if(!scalarToSolidMapping_.valid())
        {
            makeMeshMapping();
        }
        return scalarToSolidMapping_();
    }


    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    scalarDiskReader::scalarDiskReader(
        const word scalarName,
        const fvMesh &mesh,
        const fvMesh &baseMesh,
        const dictionary dict)
        : scalarName_(scalarName),
          mesh_(mesh),
          baseMesh_(baseMesh),
          dict_(dict),
          scalarToSolidMapping_(),
          FieldPtr_(),
          fieldBaseMeshPtr_(),
          fieldWasReadFromDisk_(false),
          fieldCaseDir_(dict.getOrDefault<fileName>("scalarCaseDirectory", ".")),
          fieldRunTime(Time::controlDictName,
                mesh_.time().rootPath(),
                fileName(mesh.time().caseName()/fieldCaseDir_),
                "system",
                "constant",
                true,
                true),
          fieldMesh_(IOobject
                (
                    fvMesh::defaultRegion,
                    fieldRunTime.caseConstant(),
                    fieldRunTime,
                    IOobject::MUST_READ
                )),
          curTimeIndex_(-1)
    {

        fieldRunTime.setTime(mesh_.time());
        readField();
        //Info << "reading field " << scalarName_ << " from disk" << endl;
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool Foam::scalarDiskReader::readField()
    {

            // Set fieldRunTime to be the same as the main case
            fieldRunTime.setTime(mesh_.time());

            IOobject fieldheader
            (
                scalarName_,
                fieldRunTime.timeName(),
                fieldMesh_,
                IOobject::MUST_READ
            );

            //Info << fieldheader.name()<< fieldheader.instance() << endl;

            if (fieldheader.typeHeaderOk<volScalarField>())
            {
                Info<< nl << "Reading "<< scalarName_ <<" field from time = "
                    << fieldMesh_.time().timePath() << nl << endl;

                // Read field field from field case
                const volScalarField PField(fieldheader,fieldMesh_);

                //Reset baseMesh Field so it can be read for this timestep
                fieldBaseMeshPtr_.reset();

                // Check if a single
                if (mesh_ == baseMesh_)
                {
                    FieldPtr_.clear();
                    FieldPtr_.reset
                    (
                        scalarToSolidMapping().mapSrcToTgt(PField).ptr()
                    );
                }
                else
                {
                    // See if another scalarReader from another mechanicalLaw already read
                    // the base mesh field and saved it, otherwise create it.
                    if (!baseMesh_.foundObject<volScalarField>(scalarName_))
                    {
                        // Copy PField to a field field on the baseMesh
                        fieldBaseMeshPtr_.reset(
                        new volScalarField
                        (
                            IOobject
                            (
                                scalarName_,
                                baseMesh_.time().timeName(),
                                ".", //Using region function to revert mesh region change
                                baseMesh_,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            scalarToSolidMapping().mapSrcToTgt(PField)()
                            //"zeroGradient"
                        )
                        );
                        baseMesh_.objectRegistry::checkIn(fieldBaseMeshPtr_());
                    }

                    // Map fieldBaseMesh to field on the current material mesh
                    FieldPtr_.clear();
                    FieldPtr_.reset
                    (
                        new volScalarField
                        (
                            baseMesh_.lookupObject<mechanicalModel>
                            (
                                "mechanicalProperties"
                            ).solSubMeshes().lookupBaseMeshVolField<scalar>
                            (
                                scalarName_, mesh_
                            )
                        )
                    );
                }

                fieldWasReadFromDisk_ = true;

                // Set the field field to write to disk: this allows a restart where
                // the field field was only specified in the first time-step
                // Also, it is convenient to see the field
                FieldPtr_.ref().writeOpt() = IOobject::AUTO_WRITE;
            }
            else if(mesh_.time().timeIndex()>mesh_.time().startTimeIndex())
            {
                // This time folder didnt have a field to read. So we write it.
                if(fieldBaseMeshPtr_.valid())
                {
                    Info << "writing out base-mesh field " << scalarName_ << endl;
                    fieldBaseMeshPtr_().write();
                }
                fieldWasReadFromDisk_ = false;
            }
            else
            {
                FatalError("Reading "+scalarName_) << "You have activated the scalarField reader for the field "
                                                   << scalarName_ << "in poroMaterialLaw," << nl
                                                   << "but there was no field in the start time folder." << nl
                                                   << "Did you mean to run a coupled simulation?" << endl;
            }


        return fieldWasReadFromDisk_;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
