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

\*---------------------------------------------------------------------------*/

#include "iterationControl.H"
#include "fvSolution.H"
#include "solutionControl.H"
#include <boost/format.hpp>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::iterationControl::clear()
{
    total_ = 0;
    interval_ = 0;

    onLoop_.clear();
    onConverged_.clear();
    onEnd_.clear();

    residuals_.clear();
    
    converged_ = false;
}

void Foam::iterationControl::stop()
{
    // If called manually, ensure status() will return false
    //index_ = total_ + 1;
}

void Foam::iterationControl::read(dictionary& dict)
{
    clear();

    bool enabled = dict.lookupOrAddDefault("enabled", true);

    if (enabled)
    {
        scalar timeStart;
        if (dict.readIfPresent("timeStart", timeStart))
        {
            timeStart = time_.userTimeToTime(timeStart);

            enabled =
            (
                enabled
             && time_.value() >= (timeStart - 0.5*time_.deltaTValue())
            );
        }

        scalar timeEnd;
        if (dict.readIfPresent("timeEnd", timeEnd))
        {
            timeEnd = time_.userTimeToTime(timeEnd);

            enabled =
            (
                enabled
             && time_.value() <= (timeEnd + 0.5*time_.deltaTValue())
            );
        }
    }

    if (!enabled)
    {
        return;
    }

    total_ = dict.lookupOrAddDefault("iterations", total_);
    abortCalc_ = dict.lookupOrAddDefault("maxIterAbort", abortCalc_);
    interval_ = dict.lookupOrAddDefault("interval", interval_);
    infoFrequency_ = dict.lookupOrAddDefault("infoFrequency", infoFrequency_);

    onLoop_ = dict.lookupOrAddDefault("onLoop", onLoop_);
    onConverged_ = dict.lookupOrAddDefault("onConverged", onConverged_);
    onEnd_ = dict.lookupOrAddDefault("onEnd", onEnd_);
    writeResidualFile_ = dict.lookupOrAddDefault("performanceLog", writeResidualFile_);
    Switch writeResidualField(false);
    writeResidualField = dict.lookupOrAddDefault("writeResidualField", writeResidualField);

    for (const entry& dataDictEntry : dict.subOrEmptyDict("convergence"))
    {
        residuals_.append(iterationResidual::New(time_, dataDictEntry.keyword(), dataDictEntry.stream(), writeResidualField));
    }
}

bool Foam::iterationControl::checkConverged()
{
    bool achieved = true;
    forAll(residuals_, iRes)
    {
        if( residuals_[iRes].calcResidual() < residuals_[iRes].tolerance() || residuals_[iRes].tolerance()==-1)
        {
            achieved = achieved && true;
        }
        else
        {
            achieved = false;
        }
    }
    return achieved;
}

void Foam::iterationControl::makeResidualFile()
{Info << "Creating performanceLog.dat" << endl;
    if (Pstream::master())
    {
        Info << "Creating performanceLog.dat" << endl;
        residualFilePtr_.reset(
        new OFstream(time_.path() / "performanceLog_"+name_+".dat"));
    
        residualFilePtr_()
        << "Time" << tab
        << "ExecutionTime" << tab
        << "ClockTime" << tab
        << "Iterations" << tab;
        forAll(residuals_, iRes)
        {
            residualFilePtr_() << residuals_[iRes].operationType() << residuals_[iRes].name() << tab;
        }
        residualFilePtr_() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iterationControl::iterationControl
(
    Time& runTime,
    const label nCycles,
    const word& loopName
)
:
    time_(runTime),
    index_(0),
    total_(nCycles),
    name_(loopName),
    interval_(0),
    abortCalc_(false),
    onLoop_(),
    onConverged_(),
    onEnd_(),
    infoFrequency_(1),
    converged_(false),
    residuals_(),
    residualFilePtr_(),
    writeResidualFile_(false)
{}


Foam::iterationControl::iterationControl
(
    Time& runTime,
    dictionary& algorithmDict,
    const word loopName
)
:
    iterationControl(runTime, 0, loopName)
{
        // Info<< dictName << *dictptr << endl;
        read(algorithmDict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::iterationControl::~iterationControl()
{
    residualFilePtr_.clear();
    clear();
    stop();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::iterationControl::status() const
{
    return (index_ < total_);
}

void Foam::iterationControl::reset()
{
    index_ = 0;
    converged_ = false;
    forAll(residuals_, iRes)
    {
         residuals_[iRes].reset();
    }
}


bool Foam::iterationControl::loop()
{
    bool active = (index_ < total_);   // as per status()
    if (active)
    {
        operator++();

        converged_ = checkConverged();
        info();
        if (converged_ && index_ > 1) // at least 2 iterations
        {
            time_.functionObjects().execute(onConverged_, index_);
            stop();
            return false;
        }
        else if
        (
            interval_ && !(index_ % interval_)
         && !onLoop_.empty()
        )
        {
            time_.functionObjects().execute(onLoop_, index_);
        }
    }
    else if (index_)
    {
        if(abortCalc_)
	{
	     time_.stopAt(Time::saNoWriteNow);
	}
        // Not active, the loop condition has now exiting on the last subloop
        if (!converged_ && !onEnd_.empty())
        {
            time_.functionObjects().execute(onEnd_, index_);
        }
    }
    return active;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

void Foam::iterationControl::info() const
{
    if (index() % infoFrequency_ == 0 || converged_)
                {   cout << nl << name_ << nl
                         << string(68,'=') << nl;
                    boost::format bformatr("|   %-4s %-19s  |  %15g  |   %10g   | \n");
                    cout << bformatr  %"" % "Property" % "current value" % "min/max";
                    cout << "|" << string(66,'-') << "|" << nl;
                    cout << bformatr %" " % "Time" % time_.timeOutputValue() % time_.endTime().value();
                    cout << bformatr %"#" % "Iterations" % index() % nCycles();
                    forAll(residuals_, iRes)
                    {
                        residuals_[iRes].tolerance()==-1
                        ?cout << bformatr % residuals_[iRes].operationType() % residuals_[iRes].name() % residuals_[iRes].residual() % "-"
                        :cout << bformatr % residuals_[iRes].operationType() % residuals_[iRes].name() % residuals_[iRes].residual() % residuals_[iRes].tolerance();
                    }
                    cout << string(68,'=') << nl << nl
                    << endl;
                }
}

void Foam::iterationControl::write()
{
    if(writeResidualFile_)
    {
        if(!residualFilePtr_.valid())
        {
            makeResidualFile();
        }

    residualFilePtr_.ref()
    << time_.timeName() << tab
    << time_.elapsedCpuTime() << tab
    << time_.elapsedClockTime() << tab
    << index() << tab  ;
    forAll(residuals_, iRes)
    {
        residualFilePtr_.ref() << residuals_[iRes].residual() << tab;
    }
    residualFilePtr_.ref() << endl;
    }
}

Foam::iterationControl& Foam::iterationControl::operator++()
{
    ++index_;
    return *this;
}

// ************************************************************************* //
