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
    Foam::poroSubMeshes

Description
    Class in charge of multiple sub-meshes, where each sub-mesh corresponds to
    a different material. An independent mechanical law can then be defined
    for each sub-mesh.

    Corrections are applied at bi-material interfaces to ensure continuity of
    stress.

    The sub-meshes are constructed from the cellZones of a given base mesh.

SourceFiles
    poroSubMeshes.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef poroSubMeshes_H
#define poroSubMeshes_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "newFvMeshSubset.H"
#include "newLeastSquaresVolPointInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class poroSubMeshes Declaration
\*---------------------------------------------------------------------------*/

class poroSubMeshes
{
    // Private data

        //- Const reference to the base (global) mesh
        const fvMesh& baseMesh_;

        //- Name of the cellZone for each sub-mesh
        mutable wordList cellZoneNames_;

        //- Does the solidModel use an incremental approach?
        //  i.e. does the solidModel solve for DD instead of D
        const bool incremental_;

        //- Write sub-meshes
        const bool writeSubMeshes_;

        //- Sub-meshes
        mutable PtrList<newFvMeshSubset> subMeshes_;

        //- Volume field to point field interpolator for each sub-mesh
        //mutable PtrList<newLeastSquaresVolPointInterpolation>
        //    subMeshVolToPoint_;
        /*
        //- Stress volField for each sub-mesh
        mutable PtrList<surfaceVectorField> subMeshFlux_;

        //- Stress surfaceField for each sub-mesh
        mutable PtrList<surfaceVectorField> subMeshFluxf_;

        //- Displacement (or displacement increment) volField for each
        //  sub-mesh
        mutable PtrList<volScalarField> subMeshp_;

        //- Displacement  (or displacement increment) gradient volField for each
        //  sub-mesh
        mutable PtrList<volVectorField> subMeshGradp_;

        //- Displacement  (or displacement increment) gradient surfaceField for
        //  each sub-mesh
        mutable PtrList<surfaceVectorField> subMeshGradpf_;

        //- Displacement (or displacement increment) pointField for each
        //  sub-mesh
        mutable PtrList<pointScalarField> subMeshPointp_;

        //- Flag to indicate if there are any bi-material interface faces
        mutable bool* biMaterialInterfaceActivePtr_;

        //- List of faces in the base mesh that are on a bi-material interface
        //  ie. they are shared by two materials
        mutable labelList* interfaceBaseFacesPtr_;

        //- Index of the subMesh on the other side of a bi-material interface
        mutable PtrList<labelList> interfaceShadowSubMeshID_;

        //- Index of the subMesh on the other side of a bi-material interface
        mutable PtrList<labelList> interfaceShadowPatchID_;

        //- Index of the subMesh on the other side of a bi-material interface
        mutable PtrList<labelList> interfaceShadowFaceID_;

        //- Stress calculated at shadow side of the interface for processor
        //  faces
        mutable PtrList<vectorField> interfaceShadowFlux_;

        //- Number of materials adjacent to each point in the base mesh
        mutable labelList* pointNumOfMaterialsPtr_;

        //- Isolated interface points
        mutable labelList* isolatedInterfacePointsPtr_;
        */

    // Private Member Functions

        //- Make material sub-meshes
        void makeSubMeshes() const;

        //- Make subMesh vol-to-point interpolators
//        void makeSubMeshVolToPoint() const;

        //- Check the cell zones are defined for each sub-mesh
        void checkCellZones() const;

        //- Calculate the subMesh Flux fields
//        void calcSubMeshFlux() const;

        //- Calculate the subMesh Fluxf fields
//        void calcSubMeshFluxf() const;

        //- Calculate the subMesh D fields
//        void calcSubMeshp() const;

        //- Calculate the subMesh gradD fields
//        void calcSubMeshGradp() const;

        //- Calculate the subMesh gradDf fields
//        void calcSubMeshGradpf() const;

        //- Calculate the subMesh pointD fields
//        void calcSubMeshPointp() const;

        //- Make the interface base faces
//        void makeInterfaceBaseFaces() const;

        //- Return the interface base faces
//        const labelList& interfaceBaseFaces() const;

        //- Calculate the interface shadow subMesh indices
        //  These are the maps to find which subMesh has the shadow face at a
        //  bi-material interface
//        void calcInterfaceShadowIDs() const;

        //- Return the interface shadow subMesh indices
//        const PtrList<labelList>& interfaceShadowSubMeshID() const;

        //- Return the interface shadow patch indices
//        const PtrList<labelList>& interfaceShadowPatchID() const;

        //- Return the interface shadow face indices
//        const PtrList<labelList>& interfaceShadowFaceID() const;

        //- Make materials number for points
//        void makePointNumOfMaterials() const;

        //- Return material number for points
//        const labelList& pointNumOfMaterials() const;

        //- Make isolated interface points
//        void makeIsolatedInterfacePoints() const;

        //- Return isolated interface points
//        const labelList& isolatedInterfacePoints() const;

        //- Return interfaceShadowFlux
//        const PtrList<symmTensorField>& interfaceShadowFlux() const;

        //- Make interfaceShadowFlux
//        void makeInterfaceShadowFlux() const;

        //- Assemble and update the stress calculated from the other side of the
        //  interface i.e. calculated in the subMeshes on the other side of the
        //  interface
//        void updateInterfaceShadowFlux(const bool useVolFieldFlux);

        //- Are there any bi-material interface faces
//        bool biMaterialInterfaceActive() const;

        //- Calculate biMaterialInterfaceActive
//        void calcBiMaterialInterfaceActive() const;

        //- Clear out demand driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        poroSubMeshes(const poroSubMeshes&);

        //- Disallow default bitwise assignment
        void operator=(const poroSubMeshes&);


public:

    // Constructors

        //- Construct from an fvMesh
        poroSubMeshes
        (
            const fvMesh& baseMesh,
            const wordList& cellZoneNames,
            const bool incremental,
            const bool writeSubMeshes = false
        );


    //- Destructor
    ~poroSubMeshes();


    // Member Functions

        // Access

            //- Return const access to the base mesh
            const fvMesh& baseMesh() const;

            //- Return a const reference to the subMeshes
            const PtrList<newFvMeshSubset>& subMeshes() const;

/*            //- Return a non-const reference to the subMeshes
            PtrList<newFvMeshSubset>& subMeshes();

            //- Return a const reference to the vol-to-point interpolators
            const PtrList<newLeastSquaresVolPointInterpolation>&
                subMeshVolToPoint() const;

            //- Return the subMesh Flux volFields
            PtrList<surfaceVectorField>& subMeshFlux();

            //- Return the subMesh Flux volFields
            const PtrList<surfaceVectorField>& subMeshFlux() const;

            //- Return the subMesh Flux surfaceFields
            PtrList<surfaceVectorField>& subMeshFluxf();

            //- Return the subMesh Flux surfaceFields
            const PtrList<surfaceVectorField>& subMeshFluxf() const;

            //- Return the subMesh D volFields
            PtrList<volScalarField>& subMeshp();

            //- Return the subMesh D volFields
            const PtrList<volScalarField>& subMeshp() const;

            //- Return the subMesh gradD volFields
            PtrList<volVectorField>& subMeshGradp();

            //- Return the subMesh gradD volFields
            const PtrList<volVectorField>& subMeshGradp() const;

            //- Return the subMesh gradDf volFields
            PtrList<surfaceVectorField>& subMeshGradpf();

            //- Return the subMesh gradDf volFields
            const PtrList<surfaceVectorField>& subMeshGradpf() const;

            //- Return the subMesh pointD volFields
            PtrList<pointScalarField>& subMeshPointp();

            //- Return the subMesh pointD volFields
            const PtrList<pointScalarField>& subMeshPointp() const;

            //- Lookup a field from the base mesh and interpolate it the subMesh
            template<class Type>
            tmp< GeometricField<Type, fvPatchField, volMesh> >
            lookupBaseMeshVolField
            (
                const word& fieldName,
                const fvMesh& subMesh
            ) const;

            //- Lookup a field from the base mesh and interpolate it the subMesh
            template<class Type>
            tmp< GeometricField<Type, fvsPatchField, surfaceMesh> >
            lookupBaseMeshSurfaceField
            (
                const word& fieldName,
                const fvMesh& subMesh
            ) const;

        // Edit

            //- Interpolate the base D to the subMesh D, where we apply
            //  corrections on bi-material interfaces
            void interpolateDtoSubMeshp
            (
                const volScalarField& p,
                const bool useVolFieldFlux
            );

            //- Correct snGrad on subMeshGradD fields because the default
            //  calculated patches do not account for non-orthogonality
            void correctBoundarySnGrad
            (
                PtrList<volScalarField>& subMeshDList,
                PtrList<volVectorField>& subMeshGradDList
            );

            //- Correct snGrad on subMeshGradD fields because the default
            //  calculated patches do not account for non-orthogonality
            void correctBoundarySnGradf
            (
                PtrList<volScalarField>& subMeshDList,
                PtrList<surfaceVectorField>& subMeshGradDfList,
                PtrList<volVectorField>& subMeshGradDList
            );

            //- Map a volField from the subMesh to the base mesh
            template<class Type>
            void mapSubMeshVolFields
            (
                const PtrList<GeometricField<Type, fvPatchField, volMesh> >&
                    subMeshFields,
                GeometricField<Type, fvPatchField, volMesh>& baseMeshField
            ) const;

            //- Map a surfaceField from the subMesh to the base mesh
            template<class Type>
            void mapSubMeshSurfaceFields
            (
                const PtrList
                <
                    GeometricField<Type, fvsPatchField, surfaceMesh>
                >& subMeshFields,
                GeometricField<Type, fvsPatchField, surfaceMesh>& baseMeshField
            ) const;

            //- Map pointFields from the subMesh to the base mesh
            template<class Type>
            void mapSubMeshPointFields
            (
                const PtrList
                <
                GeometricField<Type, pointPatchField, pointMesh>
                >& subMshFields,
                GeometricField<Type, pointPatchField, pointMesh>& baseMeshField
            ) const;
*/
            //- Move subMeshes
            void moveSubMeshes();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "poroSubMeshesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
