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
    poroPatchFlux

Description
    FunctionObject reports fluxes over the patch of a solid.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    poroPatchFlux.C

\*---------------------------------------------------------------------------*/

#ifndef poroPatchFlux_H
#define poroPatchFlux_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                        Class poroPatchFlux Declaration
\*---------------------------------------------------------------------------*/

    class poroPatchFlux
        : public functionObject
    {
        // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time &time_;

        //- History patch ID
        label historyPatchID_;

        //- If the history patch is found
        bool patchFound_;

        //- History file
        autoPtr<OFstream> historyFilePtr_;

        word phiName_;

        // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        poroPatchFlux(
            const poroPatchFlux &);

        //- Disallow default bitwise assignment
        void operator=(const poroPatchFlux &);

    public:
        //- Runtime type information
        TypeName("poroPatchFlux");

        // Constructors

        //- Construct from components
        poroPatchFlux(
            const word &name,
            const Time &,
            const dictionary &);

        // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop

        virtual bool execute();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary &dict);

        //- Write
        virtual bool write();
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
