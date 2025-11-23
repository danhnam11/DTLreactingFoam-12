/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Application
    foamRun

Description
    Loads and executes an OpenFOAM solver module either specified by the
    optional \c solver entry in the \c controlDict or as a command-line
    argument.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

Usage
    \b foamRun [OPTION]

      - \par -solver <name>
        Solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To run a \c rhoPimpleFoam case by specifying the solver on the
        command line:
        \verbatim
            foamRun -solver fluid
        \endverbatim

      - To update and run a \c rhoPimpleFoam case add the following entries to
        the controlDict:
        \verbatim
            application     foamRun;

            solver          fluid;
        \endverbatim
        then execute \c foamRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "solver.H"
#include "pimpleSingleRegionControl.H"
#include "setDeltaT.H"
#include "clockTime.H" //
#include "OFstream.H"  //

#include "ittnotify.h" //for profiling using vtune

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "solver",
        "name",
        "Solver name"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createMyClockTime.H" // Nam
    #include "createLogFiles.H" //Nam

    // Create VTune domain for OpenFOAM solver
    __itt_domain* domain = __itt_domain_create("OpenFOAM Solver");

    // Define named regions
    __itt_string_handle* taskTimeStep = __itt_string_handle_create("totalTimeStep");    
    __itt_string_handle* taskMomentum = __itt_string_handle_create("Step2: Momentum");
    __itt_string_handle* taskYEProp   = __itt_string_handle_create("Step3+4+5: Chem+Energy+Cal.Properties");
    __itt_string_handle* taskPressure = __itt_string_handle_create("Step6+7+8: Pressure");
    __itt_string_handle* taskUpdateProperties = __itt_string_handle_create("Step10: UpdateProperties");

    // Read the solverName from the optional solver entry in controlDict
    word solverName
    (
        runTime.controlDict().lookupOrDefault("solver", word::null)
    );

    // Optionally reset the solver name from the -solver command-line argument
    args.optionReadIfPresent("solver", solverName);

    // Check the solverName has been set
    if (solverName == word::null)
    {
        args.printUsage();

        FatalErrorIn(args.executable())
            << "solver not specified in the controlDict or on the command-line"
            << exit(FatalError);
    }
    else
    {
        // Load the solver library
        solver::load(solverName);
    }

    // Create the default single region mesh
    #include "createMesh.H"

    // Instantiate the selected solver
    autoPtr<solver> solverPtr(solver::New(solverName, mesh));
    solver& solver = solverPtr();

    // Create the outer PIMPLE loop and control structure
    pimpleSingleRegionControl pimple(solver.pimple);

    // Set the initial time-step
    setDeltaT(runTime, solver);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        __itt_task_begin(domain, __itt_null, __itt_null, taskTimeStep);
            
        #include "createTimers.H" //Nam
        // Update PIMPLE outer-loop parameters if changed
        Time_pre = myTime.elapsedTime(); // timer 0
        pimple.read();

        solver.preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solver);

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;
        Time0 = myTime.elapsedTime(); // timer 0

        // PIMPLE corrector loop
        while (pimple.loop())
        {
            solver.moveMesh();
            solver.motionCorrector();
            solver.fvModels().correct();

            Time1 = myTime.elapsedTime(); // timer 1
            solver.prePredictor();

            __itt_task_begin(domain, __itt_null, __itt_null, taskMomentum);
            Time2 = myTime.elapsedTime(); // timer 2
            solver.momentumPredictor();
            __itt_task_end(domain);

            __itt_task_begin(domain, __itt_null, __itt_null, taskYEProp);
            Time3 = myTime.elapsedTime(); // timer 3           
            solver.thermophysicalPredictor();
            __itt_task_end(domain);            

            __itt_task_begin(domain, __itt_null, __itt_null, taskPressure);
            Time4 = myTime.elapsedTime(); // timer 4       
            solver.pressureCorrector();
            __itt_task_end(domain);

            Time5 = myTime.elapsedTime(); // timer 5       
            solver.postCorrector();
            //Time6 = myTime.elapsedTime(); // timer 6       
        }

        __itt_task_begin(domain, __itt_null, __itt_null, taskUpdateProperties);        
        Time7 = myTime.elapsedTime(); // timer 7       
        //- Update properties for nex time step - Nam 
        solver.postSolve();
        __itt_task_end(domain);

        Time8 = myTime.elapsedTime(); // timer 8
        #include "calculateTimeAndWriteLogFiles.H" // Nam

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        __itt_task_end(domain);
            
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
