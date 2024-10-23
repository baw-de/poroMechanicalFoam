

DSigma_(
    IOobject(
        "DSigma",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedSymmTensor("DSigma", dimForce / dimArea, symmTensor::zero)),
    P_(
        IOobject(
            "P",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("P", dimPressure, 0.0)),
    s_(
        IOobject(
            "s",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedSymmTensor("s", dimForce / dimArea, symmTensor::zero)),

    r_(
        IOobject(
            "r",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("r", dimless, symmTensor::zero)),

    alpha_(
        IOobject(
            "alpha",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("alpha", dimless, symmTensor::zero)),
    F_(
        IOobject(
            "F",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("F", dimless, symmTensor::zero)),
    m_(
        IOobject(
            "m",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("m", dimless, m0_.value())),

    DEpsilon_(
        IOobject(
            "DEpsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("DEpsilon", dimless, symmTensor::zero)),

    DEpsilonV_(
        IOobject(
            "DEpsilonV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("DEpsilonV", dimless, 0.0)),
    DEpsilonD_(
        IOobject(
            "DEpsilonD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("DEpsilonD", dimless, symmTensor::zero)),
    epsilonV_(
        IOobject(
            "epsilonV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("epsilonV", dimless, 0.0)),
    epsilonD_(
        IOobject(
            "epsilonD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedSymmTensor("epsilonD", dimless, symmTensor::zero)),

    DEpsilonP_(
        IOobject(
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)),
    DEpsilonPV_(
        IOobject(
            "DEpsilonPV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("DEpsilonPV", dimless, 0.0)),
    DEpsilonPD_(
        IOobject(
            "DEpsilonPD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor("DEpsilonPD", dimless, symmTensor::zero)),
    epsilonP_(
        IOobject(
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)),
    e_(
        IOobject(
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("e", dimless, e0_.value())),

    psi_(
        IOobject(
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("psi", dimless, (e0_ - e_cref_).value())),
    K_(
        IOobject(
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("K", dimPressure, K0_.value())),

    mu_(
        IOobject(
            "G",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("G", dimPressure, G0_.value())),
