autoPtr<meshToMesh> solidToPoroFluid_;

    // Attempt to read the poroCouplingProperties dictionary
    IOdictionary meshToMeshDict
    (
        IOobject
        (
            "poroCouplingProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Get the mapping method
    const word mapMethodName
    (
        meshToMeshDict.lookupOrDefault<word>("mapMethod", "direct")
    );

    // Check mapping method is valid
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << "Unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalError);
    }

    // Convert string to enum value
    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    // Check if consistent mapping is required
    Switch consistent
    (
        meshToMeshDict.lookupOrDefault<Switch>("consistent", true)
    );

    if (!consistent)
    {
        if (mapMethodName == "direct")
        {
            WarningInFunction
                << "Direct mapping selected for inconsistent region meshes, "
                << "changing to imMapNearest"
                << endl;
            mapMethod = meshToMesh::interpolationMethod::imMapNearest;
        }

        // Read patch mapping settings
        HashTable<word> patchMap;
        wordList cuttingPatches;

        meshToMeshDict.readEntry("patchMap", patchMap);
        meshToMeshDict.readEntry("cuttingPatches", cuttingPatches);

        // Create the mapping object
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
        // Create the mapping object with default settings
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
