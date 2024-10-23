/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::iterationControl

Description
    A class for managing arbitrary loops, calculate and check convergence criteria, generate
    a convergence metric output file and execute function objects.

Usage
    Examples of function object specification:
    \verbatim
    SIMPLE
    {
        energyCoupling
        {
            iterations  100;
            onLoop      ();
            onConverged ( externalCoupled  "loopThings.*" );

            convergence
            {
                // Option 1: L2 norm of mag(D-D.prevIter())
                "D"  L2 1e-3;
                // Option 2: maximum of D-D.prevIter()/max(mag(D))
                "rel D" max 1e-3;
                // Option 3: Calculate and show root-mean-squared residual, don't check for convergence
                "D" RMS show;
                // Option 4: Use linearSolver residuals
                linearSolver
            }
        }
    }
    \endverbatim

    Where the loop entries comprise:
    \table
        Property    | Description           | Required      | Default
        enabled     | active/deactive loop  | no            | true
        iteration   | times to loop         | no            | 0
        timeStart   | begin time for loop activation  | no  | -VGREAT
        timeEnd     | end time of loop activation     | no  | VGREAT
        interval    | sub-interval to execute onLoop  | no  | 0
        onLoop      | function object names to call at executeInterval | no
        onConverged | function object names to call when converged | no
        onEnd       | function object names to call when loop ends | no
        convergence | dictionary of convergence values to check | no
    \endtable

    The function object names listed by \c onLoop, \c onConverged, \c onEnd
    must implement an \c execute(int) method.
    If the time controls \c timeStart or \c timeEnd are used for the loop,
    these values are only inspected upon creation, not during execution.
    
Author
    Denis Maier, BAW. All rights reserved.


SeeAlso
    fvSolution

Author
    Denis Maier, BAW. All Rights Reserved.

SourceFiles
    iterationControl.C

\*---------------------------------------------------------------------------*/

#ifndef iterationControl_H
#define iterationControl_H

#include "subLoopTime.H"
#include "iterationResidual.H"
#include "dictionary.H"
#include "wordRes.H"
#include "OFstream.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class iterationControl Declaration
\*---------------------------------------------------------------------------*/

class iterationControl
{
    // Private Member Functions

        //- Clear
        void clear();

        void stop();

        //- Read settings from dictionary
        void read(dictionary& dict);

        bool checkConverged();

        void makeResidualFile();

        //- No copy construct
        iterationControl(const iterationControl&) = delete;

        //- No copy assignment
        void operator=(const iterationControl&) = delete;

protected:

    // Protected data

        //- reference to the time object
        Time& time_;

        //- The current sub-cycle index
        label index_;

        //- The total number of cycles in the sub-cycle
        label total_;

        //- Name of the loop control (the lookup dictionary name).
        word name_;

        //- The interval to execute onLoop function-objects
        label interval_;

	    //- Flag what to do when max. iterations are reached
	    Switch abortCalc_;

        //- Function object names to fire during the loop (at executeInterval)
        wordRes onLoop_;

        //- Function object names to fire on convergence
        wordRes onConverged_;

        //- Function object names to fire when the loop exits without
        //- convergence
        wordRes onEnd_;

        //- Frequency of residual information output
        label infoFrequency_;

        //- Convergence tests passed
        bool converged_;

        //- List of residuals to be compared
        PtrList<iterationResidual> residuals_;

        // residual output logfile
        mutable autoPtr<OFstream> residualFilePtr_;

        //- write Residual file
        Switch writeResidualFile_;

        //- Switch to write the residual field
        Switch writeResidualField_;

public:

    // Constructors

        //- Construct from time with fixed number of cycles
        //  \param runTime  the top-level time
        //  \param nCycles  the number of times to loop
        //  \param loopName  the name of the loop
        iterationControl
        (
            Time& runTime,
            const label nCycles,
            const word& dictName = "loop"
        );

        //- Construct from fvSolution dictionary based on time and the name
        //- of the controlling algorithm
        //  \param runTime  the top-level time
        //  \param algorithmDict the fvSolution algorithm dictionary,
        //      typically PIMPLE or SIMPLE
        //  \param dictName  the name of the control dictionary
        iterationControl
        (
            Time& runTime,
            dictionary& algorithmDict,
            const word loopName = "loop"
        );


    //- Destructor
    ~iterationControl();


    // Member Functions

        //- Name of the loop control
        inline const word& name() const
        {
            return name_;
        }

        //- The current cycle index
        inline label index() const
        {
            return index_;
        }

        //- The total number of cycles
        inline label nCycles() const
        {
            return total_;
        }

        //- The interval to execute onLoop function-objects
        inline label interval() const
        {
            return interval_;
        }

        PtrList<iterationResidual>& residuals()
        {
            if(residuals_.size()==0)
            {
                WarningInFunction << "No residuals are selected!" << endl;
            }
            return residuals_;
        }

        //- Return true if looping is active
        bool status() const;

        //- reset the loop
        void reset();

        //- True if looping is active, increments the index and executes
        //- the onLoop and onConverged functions.
        //  Example usage,
        //  \code
        //      while (control.loop())
        //      {
        //          solve;
        //      }
        //  \endcode
        bool loop();

        void info() const;

        void write();

        // Member Operators

        //- Prefix increment
        iterationControl& operator++();
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
