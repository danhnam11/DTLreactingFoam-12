/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "profilingDTLreactingFoam.H"
#include "fvcDdt.H"

#include "ittnotify.h" //for profiling using vtune

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::profilingDTLreactingFoam::thermophysicalPredictor()
{

    // Create VTune domain for thermophysicalPredictor function
    __itt_domain* domain = __itt_domain_create("inside Step345");

    // Define named regions
    __itt_string_handle* taskChem = __itt_string_handle_create("Step3: Chemistry");
    __itt_string_handle* taskEnergy = __itt_string_handle_create("Step4: Energy");
    __itt_string_handle* taskProperties = __itt_string_handle_create("Step5: Cal.Properties");

    // create myTime object for precise time recording
    Foam::clockTime myTime;

    // create timers
    scalar Time0, Time1,Time2,Time3,Time4,Time5, Time6;
    Time0=0; Time1 = 0;Time2 = 0;Time3 = 0;Time4 = 0; Time5 = 0; Time6=0;

    Time0 = myTime.elapsedTime(); // timer 0
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.schemes().div("div(phi,Yi_h)")
        )
    );

    __itt_task_begin(domain, __itt_null, __itt_null, taskChem);

    Time1 = myTime.elapsedTime(); // timer 1
    reaction->correct();
    Time2 = myTime.elapsedTime(); // timer 2


    forAll(Y, i)
    {
        volScalarField& Yi = Y_[i];

        // calculate correction terms - Nam 
        sumYDiffusionCorrections1 *= 0.;
        sumYDiffusionCorrections2 *= 0.;
        forAll(Y, k)
        {
            sumYDiffusionCorrections1 += fvc::div(Y[i]*rho*Dimix[k]* fvc::grad(Y[k]));
            sumYDiffusionCorrections2 += fvc::div(Y[i]*rho*Dimix[k]*Y[k]/Wmix * fvc::grad(Wmix));
        }
        //

        if (thermo_.solveSpecie(i))
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              //+ thermophysicalTransport->divj(Yi)
              - fvm::laplacian(rho*Dimix[i], Yi) //
             ==
                reaction->R(Yi)
              //+ fvc::laplacian(rho*Dimix[i]*Yi/Wmix, Wmix) //
              + fvc::div(rho*Dimix[i]*Yi/Wmix * fvc::grad(Wmix)) //              
              - sumYDiffusionCorrections1 //
              - sumYDiffusionCorrections2 //                     
              + fvModels().source(rho, Yi)
            );

            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
        else
        {
            Yi.correctBoundaryConditions();
        }
    }

    thermo_.normaliseY();

    Time3 = myTime.elapsedTime(); // timer 3

    // Nam 
    //-Calculate diffusion velocity of inert species
    volVectorField YVt(0.0*YVi[0]);

    forAll(YVi, i)
    {
        //- Update diffusion velocity for all species
        YVi[i] = (-(Dimix[i]/Wmix)*Wmix*fvc::grad(Y[i])-(Dimix[i]/Wmix)*Y[i]*fvc::grad(Wmix));
        if (thermo_.solveSpecie(i))
        {
            YVt += YVi[i];
        }
    }
    YVi[thermo_.defaultSpecie()] = -YVt;  
    // Nam - end

    Time4 = myTime.elapsedTime(); // timer 4

    __itt_task_end(domain);    


    __itt_task_begin(domain, __itt_null, __itt_null, taskEnergy);

    volScalarField& he = thermo_.he();
    volScalarField& T = thermo_.T(); // Nam

    // Nam 
    //-Reset and then calculate Heat by diffusion terms-for each time step
    sumHeatDiffusion1 *= 0.;
    sumHeatDiffusion2 *= 0.;

    forAll(Y, k)
    {
        sumHeatDiffusion1 += fvc::div((kappa/Cp)*hei[k] * fvc::grad(Y[k]));
        sumHeatDiffusion2 += fvc::div(hei[k]*rho*YVi[k]);
    }

    heatOfReactions = reaction->Qdot();
    // Nam 

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + pressureWork
        (
            he.name() == "e"
          ? mvConvection->fvcDiv(phi, p/rho)()
          : -dpdt
        )
      //+ thermophysicalTransport->divq(he)
      - fvm::laplacian(kappa/Cp, he) //      
     ==
        reaction->Qdot()
      - sumHeatDiffusion1 //
      - sumHeatDiffusion2 //   
      + (
            buoyancy.valid()
          ? fvModels().source(rho, he) + rho*(U & buoyancy->g)
          : fvModels().source(rho, he)
        )
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(he);

    Time5 = myTime.elapsedTime(); // timer 5

    __itt_task_end(domain);    

    __itt_task_begin(domain, __itt_null, __itt_null, taskProperties);
    // for coTHERM
    calculateAndUpdateCoTHERMFlags(); 
    thermo_.correct();

    Time6 = myTime.elapsedTime(); // timer 6

    __itt_task_end(domain);    


    // Nam, monitor the temperature for more convenience
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;



    // calculate CPU times for single steps
    scalar chemTime = Time2 - Time1; 
    scalar YiTime = Time3 - Time2; 
    scalar VYiTime = Time4 - Time3;
    scalar ETime = Time5 - Time4; 
    scalar calPropertiesTime = Time6 - Time5; 
    scalar totalYEcalPropertiesTime =Time6 - Time0;
    scalar timeStep3 = chemTime + YiTime + VYiTime;

    // write the calculated CPU times to logFile 
    logFile 
        << runTime.value() << ","
        << chemTime << ","
        << YiTime << ","
        << VYiTime << ","
        << ETime << ","
        << calPropertiesTime << ","
        << totalYEcalPropertiesTime << endl;


// show CPU time measurement on the screen
Info << "CPU time inside Step345 of Time = " << runTime.value() << endl;
Info << tab << "Step3 = " << timeStep3 << " (s)," 
    << tab << 100*timeStep3/totalYEcalPropertiesTime << "[%], Chemistry" << endl;
Info << tab << "Step4 = " << ETime << " (s)," 
    << tab << 100*ETime/totalYEcalPropertiesTime << "[%], Energy" << endl;
Info << tab << "Step5 = " << calPropertiesTime << " (s)," 
    << tab << 100*calPropertiesTime/totalYEcalPropertiesTime << "[%], Cal.Properties" << endl;
Info << tab << "total = " << totalYEcalPropertiesTime << " (s)," 
    << tab << 100*totalYEcalPropertiesTime/totalYEcalPropertiesTime << "[%], totalTime_step345" << endl;

}


// ************************************************************************* //
