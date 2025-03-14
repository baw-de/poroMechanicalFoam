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

\*----------------------------------------------------------------------------*/

#include "poroPatchFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroPatchFlux, 0);

    addToRunTimeSelectionTable(
        functionObject,
        poroPatchFlux,
        dictionary);
} // namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::poroPatchFlux::writeData()
{
    if (patchFound_)
    {
        // Lookup the solid mesh
        const fvMesh *meshPtr = NULL;
        if (time_.foundObject<fvMesh>("solid"))
        {
            meshPtr = &(time_.lookupObject<fvMesh>("solid"));
        }
        else
        {
            meshPtr = &(time_.lookupObject<fvMesh>("region0"));
        }
        const fvMesh &mesh = *meshPtr;

        // Patch area vectors
        /*const vectorField& patchSf =
            mesh.Sf().boundaryField()[historyPatchID_];*/

        // Patch unit area vectors
        const vectorField patchNf = mesh.boundary()[historyPatchID_].nf();

        // Calculate the flux as the intergal over the area
        scalar flux = 0.0;

        // Lookup the stress field
        const scalarField &phi =
            mesh.lookupObject<surfaceScalarField>(
                    phiName_)
                .boundaryField()[historyPatchID_];

        // Linear geometry

        flux = gSum(phi);

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                << " " << flux
                << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poroPatchFlux::poroPatchFlux(
    const word &name,
    const Time &t,
    const dictionary &dict)
    : functionObject(name),
      name_(name),
      time_(t),
      historyPatchID_(-1),
      patchFound_(false),
      historyFilePtr_(),
      phiName_(dict.lookupOrDefault<word>("flux", "phi"))
{
    Info << "Creating " << this->name() << " function object" << endl;

    word historyPatchName("notSpecified");

    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch") >> historyPatchName;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "poroPatchFlux: historyPatch not specified" << endl;
    }

    // Lookup the solid mesh
    const fvMesh *meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh &mesh = *meshPtr;

    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);

    if (historyPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found"
            << endl;
    }
    else
    {
        patchFound_ = true;
    }

    // Create history file if not already created
    if (!historyFilePtr_.valid() && patchFound_)
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path() / ".." / "history" / startTimeName;
            }
            else
            {
                historyDir = time_.path() / "history" / startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset(
                new OFstream(
                    historyDir / "poroPatchFlux" + historyPatchName + ".dat"));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time flux" << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::poroPatchFlux::start()
{
    return writeData();
}

bool Foam::poroPatchFlux::execute()
{
    return writeData();
}

bool Foam::poroPatchFlux::read(const dictionary &dict)
{
    return true;
}

bool Foam::poroPatchFlux::write()
{
    return writeData();
}

// ************************************************************************* //
