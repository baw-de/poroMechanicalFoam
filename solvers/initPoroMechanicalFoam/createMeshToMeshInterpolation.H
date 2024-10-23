    autoPtr<meshToMesh> solidToPoroFluid_;

    IOdictionary meshToMeshDict(
          IOobject(
              "poroCouplingProperties",
              runTime.constant(),
              runTime,
              IOobject::MUST_READ,
              IOobject::NO_WRITE));

    const word mapMethodName(meshToMeshDict.lookupOrAddDefault<word>("mapMethod","direct"));
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
          meshToMeshDict.lookupOrAddDefault<Switch>(
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

        meshToMeshDict.readEntry("patchMap", patchMap);
        meshToMeshDict.readEntry("cuttingPatches", cuttingPatches);
        solidToPoroFluid_.reset
                (
                    new meshToMesh
                    (
                        fluid().mesh(),
                        solid().mesh(),
                        mapMethod,
                        patchMap,
                        cuttingPatches
                    )
                );
        }
        else
        {
        solidToPoroFluid_.reset
                (
                    new meshToMesh
                    (
                        fluid().mesh(),
                        solid().mesh(),
                        mapMethod
                    )
                );
        }
