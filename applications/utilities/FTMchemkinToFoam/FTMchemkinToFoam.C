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
    FTMchemkinToFoam

    developed by Danh Nam Nguyen and Jae Hun Lee,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description

    Utility to generate input file for using DTLreactingFoam with FTM 

Usage

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "solver.H"
#include "pimpleSingleRegionControl.H"
#include "clockTime.H"
#include "OFstream.H"
#include "fluidMulticomponentThermo.H"
#include "psiMulticomponentThermo.H"
#include "DTMMulticomponentMixture.H"
#include "specie.H"
#include "janafThermo.H"
#include "perfectGas.H"
#include "sensibleEnthalpy.H"
#include "DTMTransport.H"
#include "preprocessingFTMTransport.H"
#include "psiThermo.H"
#include "PsiThermo.H"
#include "thermo.H"
#include "BasicThermo.H"
#include "basicThermo.H"
#include "GeometricFieldListSlicer.H"

#include "DPOLFT.H"
#include "DPCOEF.H"
#include "DP1VLU.H"
#include "PropertyReader.H"

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
    // Create the default single region mesh
    #include "createMesh.H"
    #include "createFields.H"
    #include "FittingFTM.H"
    #include "WriteFTM.H"

    Info << "Patched thermo.DTM → thermo.FTM" << endl;

    Info<< "Dij of N2 = " << mixture_.Dij(0,1, 10325, 300) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
