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

Class
    Foam::solvers::profilingDTLreactingFoam

Description
    Solver module for steady or transient turbulent flow of compressible
    multicomponent fluids with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, chemical reactions,
    combustion, Lagrangian particles, radiation, surface film etc. and
    constraining or limiting the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    profilingDTLreactingFoam.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef profilingDTLreactingFoam_H
#define profilingDTLreactingFoam_H

#include "isothermalFluid.H"
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"

#include "clockTime.H" //
#include "OFstream.H"  //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class profilingDTLreactingFoam Declaration
\*---------------------------------------------------------------------------*/

class profilingDTLreactingFoam
:
    public isothermalFluid
{

protected:

    // Thermophysical properties
        fluidMulticomponentThermo& thermo_;

    // Composition
        PtrList<volScalarField>& Y_;

    // Reactions
        autoPtr<combustionModel> reaction;
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    // Thermophysical transport
        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport;

    // Fields for DTM 
        // Molecular weight of mixture, [kg/kmol] 
        volScalarField Wmix;
        // Dynamic viscosity of mixture [kg/m.s], for laminar 
        volScalarField mu;
        // Constant pressure heat capacity, [J/kgK] 
        volScalarField Cp;
        // Thermo conductivity of mixture [W/mK], for laminar 
        volScalarField kappa;
        // Enthalpy/Internal field of mixture [J/kg]
        volScalarField HE;
        // Correction terms for diffusion velocity
        volScalarField sumYDiffusionCorrections1;
        volScalarField sumYDiffusionCorrections2;
        // Terms in energy equation         
        volScalarField sumHeatDiffusion1;
        volScalarField sumHeatDiffusion2;
        // Heat of reaction 
        volScalarField heatOfReactions;
        // List of difusion velocity and mass fraction of species
        PtrList<volVectorField> YVi;
        //- List of mass diffusivity fields [m^2/s] of species in the mixture
        PtrList<volScalarField> Dimix;
        //- List of individual Enthalpy/Internal Energy fields [J/kg] of species        
        PtrList<volScalarField> hei;

private:

    // Private Member data
        IOdictionary thermDict;
        List<word> CoTHERMMajorSpecies;
        scalar epsilonT;        
        scalar epsilonSpecies;
        List<label> CoTHERMSpeciesIndex;
        volScalarField flagSpecies;

        OFstream logFile;

    // Private Member Functions
        //- Set rDeltaT for LTS
        virtual void setRDeltaT();

        //- Calculate and update flags for coTHERM
        void calculateAndUpdateCoTHERMFlags();

public:

    // Public Data
        //- Reference to the fluid thermophysical properties
        const fluidMulticomponentThermo& thermo;

        //- Reference to the composition
        const PtrList<volScalarField>& Y;

        //- we don't need public access to the fields for DTM - Nam


    //- Runtime type information
    TypeName("profilingDTLreactingFoam");


    // Constructors

        //- Construct from region mesh
        profilingDTLreactingFoam(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        profilingDTLreactingFoam(const profilingDTLreactingFoam&) = delete;


    //- Destructor
    virtual ~profilingDTLreactingFoam();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step - Nam 
        virtual void postSolve();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const profilingDTLreactingFoam&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
