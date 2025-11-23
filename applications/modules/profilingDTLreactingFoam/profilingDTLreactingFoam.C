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

\*---------------------------------------------------------------------------*/

#include "profilingDTLreactingFoam.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"
#include <vector> //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(profilingDTLreactingFoam, 0);
    addToRunTimeSelectionTable(solver, profilingDTLreactingFoam, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::profilingDTLreactingFoam::profilingDTLreactingFoam(fvMesh& mesh)
:
    isothermalFluid
    (
        mesh,
        autoPtr<fluidThermo>(fluidMulticomponentThermo::New(mesh).ptr())
    ),

    thermo_(refCast<fluidMulticomponentThermo>(isothermalFluid::thermo_)),

    Y_(thermo_.Y()),

    reaction(combustionModel::New(thermo_, momentumTransport())),

    thermophysicalTransport
    (
        fluidMulticomponentThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo_
        )
    ),
    // Nam - for DTM
    Wmix
    (
        IOobject
        (
            "Wmix",
            runTime.name(),
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE 
        ),
        thermo_.Wmix()
    ),
    mu
    (
        IOobject
        (
            "mu",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.mu()
    ),
    Cp
    (
        IOobject
        (
            "Cp",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.Cp()
    ),
    kappa
    (
        IOobject
        (
            "kappa",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.kappa()
    ),
    HE
    (
        IOobject
        (
            "HE",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),  
        thermo_.he()
    ),
    sumYDiffusionCorrections1
    (
        IOobject
        (
            "gas_sumCorrections1",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("gas_sumCorrections1", dimensionSet(1, -3, -1, 0, 0), 0.)
    ),
    sumYDiffusionCorrections2
    (
        IOobject
        (
            "gas_sumCorrections2",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("gas_sumCorrections2", dimensionSet(1, -3, -1, 0, 0), 0.)
    ),
    sumHeatDiffusion1
    (
        IOobject
        (
            "sumHeatDiffusion",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sumHeatDiffusion", dimensionSet(1, -1, -3, 0, 0), 0.)
    ),
    sumHeatDiffusion2
    (
        IOobject
        (
            "sumHeatDiffusion2",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sumHeatDiffusion2", dimensionSet(1, -1, -3, 0, 0), 0.)
    ), 
    heatOfReactions
    (
        IOobject
        (
            "heatOfReactions",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("heatOfReactions", dimensionSet(1, -1, -3, 0, 0), 0.)
    ), 
    YVi(Y_.size()),
    Dimix(Y_.size()),
    hei(Y_.size()),
    // Nam - for coTHERM
    thermDict
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    CoTHERMMajorSpecies(thermDict.lookup("majorSpeciesForCoTHERM")),
    epsilonT(thermDict.lookupOrDefault("epsilonT", 0.1)),
    epsilonSpecies(thermDict.lookupOrDefault("epsilonS", 0.001)),
    CoTHERMSpeciesIndex(CoTHERMMajorSpecies.size()),
    flagSpecies
    (
        IOobject
        (
            "flagSpecies",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("flagSpecies", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    logFile
    (
        runTime.path()/("logCPUTime_YEcalProperties_" + runTime.userTimeName() + ".csv")
    ),
    // Nam
    thermo(thermo_),
    Y(Y_)
{
    thermo.validate(type(), "h", "e");

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());

    // Nam - for DTM
    forAll(YVi, i) 
    {  
        YVi.set
        (
            i, 
            new volVectorField 
            ( 
                IOobject
                (
                "YVi_"+thermo_.species()[i],
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE 
                ), 
                mesh,
                dimensionedVector("YVi_"+thermo_.species()[i], dimensionSet(0,1,-1, 0, 0), vector(0.,0.,0.))
            )
        ); 
    }

    forAll(Dimix, i) 
    {  
        Dimix.set
        (
            i, 
            new volScalarField 
            ( 
                IOobject
                (
                "Dimix_"+thermo_.species()[i],
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE 
                ), 
                thermo_.Dimix(i)
            )
        ); 
    }

    forAll(hei, i) 
    {  
        hei.set
        (
            i, 
            new volScalarField 
            ( 
                IOobject
                (
                "hei_"+thermo_.species()[i],
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE 
                ), 
                thermo_.hei(i)
            )
        ); 
    }

    // Nam -for coTHERM
    forAll(CoTHERMSpeciesIndex, i)
    {
        if (!thermo_.species().found(CoTHERMMajorSpecies[i]))
        {
            Info << "\n[Error***!], cannot find " << CoTHERMMajorSpecies[i] << " in the chemical Mechanism \n" << endl;
            Info << "Please, adjust the 'majorSpeciesForCoTHERM' in 'physicalProperties' dictionary!"  << exit(FatalIOError);
        }
        else
        {
            CoTHERMSpeciesIndex[i] = thermo_.species()[CoTHERMMajorSpecies[i]];
        }
    }
    Info << "\nMajor species are selected for CoTHERM are " << CoTHERMMajorSpecies << endl;
    Info << "epsilon for Temperature in CoTHERM = " << epsilonT << " [K]" << endl;
    Info << "epsilon for Species in CoTHERM     = " << epsilonSpecies << " [-](in mass fraction)\n" << endl;

    logFile
        << "Time" << ","  
        << "chemTime(s)" << ","
        << "YiTime(s)" << ","
        << "VYiTime(s)" << ","
        << "ETime(s)" << ","
        << "calPropertiesTime(s)" << ","
        << "totalYEcalPropertiesTime(s)" << endl;   
    // Nam
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::profilingDTLreactingFoam::~profilingDTLreactingFoam()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::profilingDTLreactingFoam::prePredictor()
{
    isothermalFluid::prePredictor();

    if (pimple.predictTransport())
    {
        thermophysicalTransport->predict();
    }
}


void Foam::solvers::profilingDTLreactingFoam::postCorrector()
{
    isothermalFluid::postCorrector();

    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
}


// New postSolver function for DTM
void Foam::solvers::profilingDTLreactingFoam::postSolve()
{
    //- Update properties for nex time step
    // Info << "This is postSolver() function in profilingDTLreactingFoam class" << endl;
    Info << "Update properties for next time step or iteration" << endl;

    isothermalFluid::postSolve();

    mu = thermo.mu();
    Cp = thermo.Cp();
    kappa = thermo.kappa();
    Wmix = thermo.Wmix();
    HE = thermo.he();
    forAll(Dimix, i)
    {
         Dimix[i] = thermo.Dimix(i);
         hei[i] = thermo.hei(i);
    }
}


//- Calculate and update flags for coTHERM
void Foam::solvers::profilingDTLreactingFoam::calculateAndUpdateCoTHERMFlags()
{
    label nSCoTHERM(CoTHERMMajorSpecies.size()); 
    
    //Info << "\nCalculating flags for coTHERM" << endl; 
    // calculate flags for internal fields
    forAll(flagSpecies, celli)
    {
        //initialize the flag for species in a single cell 
        std::vector<double> flagSpeciesVector(nSCoTHERM, 1.0);  
        forAll(CoTHERMMajorSpecies, j)
        {
    
            if 
            (
                mag
                (
                    Y_[CoTHERMSpeciesIndex[j]].oldTime()[celli] 
                - Y_[CoTHERMSpeciesIndex[j]][celli]
                ) <= epsilonSpecies
            )
            {
                flagSpeciesVector[j] = 0.0;
            }
            else
            {
                flagSpeciesVector[j] = 1.0;
            }
        }
    
        if 
        (
            std::all_of
            (
                flagSpeciesVector.begin(), 
                flagSpeciesVector.end(), 
                [&](double v){ return v == 0.0; }
            )
        )
        {
            // Info << "flags of all species are equal " << endl;
            flagSpecies[celli] = 0.0;
        }
        else
        {
            // Info << "flags of all species are not equal:-> need to re-calculate Properties " << endl;    
            flagSpecies[celli] = 1.0;
        }
    }
    
    // calculate flags for boundary fields
    volScalarField::Boundary& pFlagSpecies = flagSpecies.boundaryFieldRef();
    
    forAll(pFlagSpecies, patchi)
    {
        forAll(pFlagSpecies[patchi], facei)
        {
            //initialize the flag for species in a single face 
            std::vector<double> flagSpeciesVectorF(nSCoTHERM, 1.0);
            forAll(CoTHERMMajorSpecies, j)
            {
                if 
                ( 
                    mag
                    (
                        Y_[CoTHERMSpeciesIndex[j]].oldTime().boundaryField()[patchi][facei] 
                    - Y_[CoTHERMSpeciesIndex[j]].boundaryField()[patchi][facei]
                    ) <= epsilonSpecies
                )
                {
                flagSpeciesVectorF[j] = 0.0;
                }
                else
                {
                flagSpeciesVectorF[j] = 1.0;
                }
            }
    
            if 
            (
                std::all_of
                (
                    flagSpeciesVectorF.begin(), 
                    flagSpeciesVectorF.end(), 
                    [&](double v){ return v == 0.0; }
                )
            )
            {
                // Info << "flags of all species are equal " << endl;
                pFlagSpecies[patchi][facei] = 0.0;
            }
            else
            {
                // Info << "flags of all species are not equal: -> need to re-calculate Properties " << endl;   
                pFlagSpecies[patchi][facei] = 1.0;
            }
        }
    }
    
    // update flags 
    thermo_.updateCoTHERMFlags(flagSpecies); 
    //Info << "End of calculating flags for coTHERM\n" << endl; 
}


// ************************************************************************* //
