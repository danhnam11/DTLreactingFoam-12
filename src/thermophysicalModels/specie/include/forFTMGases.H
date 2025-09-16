/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#ifndef forFTMGases_H
#define forFTMGases_H

#include "specie.H"

#include "perfectGas.H"

#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "FTMTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forFTMGasEqns(Mu, He, Cp, Macro, Args...)                              \
    forThermo(Mu, He, Cp, perfectGas, specie, Macro, Args)

#define forFTMGasEnergiesAndThermos(Mu, Macro, Args...)                        \
    forFTMGasEqns(Mu, sensibleEnthalpy, janafThermo, Macro, Args);             \
    forFTMGasEqns(Mu, sensibleInternalEnergy, janafThermo, Macro, Args)

#define forFTMGasTransports(Macro, Args...)                                    \
    forFTMGasEnergiesAndThermos(FTMTransport, Macro, Args)

#define forFTMGases(Macro, Args...)                                            \
    forFTMGasTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
