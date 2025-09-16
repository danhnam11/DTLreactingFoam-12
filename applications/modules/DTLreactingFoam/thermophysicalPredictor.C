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

#include "DTLreactingFoam.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DTLreactingFoam::thermophysicalPredictor()
{
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

    reaction->correct();

    forAll(Y, i)
    {
        volScalarField& Yi = Y_[i];

        // calculate correction terms
        sumYDiffusionCorrections1 *= 0.;
        sumYDiffusionCorrections2 *= 0.;
        forAll(Y, k)
        {
            sumYDiffusionCorrections1 += fvc::div(Y[i]*rho*Dimix[k]* fvc::grad(Y[k]));
            sumYDiffusionCorrections2 += fvc::div(Y[i]*rho*Dimix[k]*Y[k]/Wmix * fvc::grad(Wmix));
        }

        if (thermo_.solveSpecie(i))
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(rho*Dimix[i], Yi)
             ==
                reaction->R(Yi)
              + fvc::div(rho*Dimix[i]*Yi/Wmix * fvc::grad(Wmix))             
              - sumYDiffusionCorrections1
              - sumYDiffusionCorrections2                  
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

    volScalarField& he = thermo_.he();
    volScalarField& T = thermo_.T();

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
      - fvm::laplacian(kappa/Cp, he)    
     ==
        reaction->Qdot()
      - sumHeatDiffusion1
      - sumHeatDiffusion2
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

    // for coTHERM
    calculateAndUpdateCoTHERMFlags(); 
    
    thermo_.correct();

    // monitor the temperature for more convenience
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}


// ************************************************************************* //
