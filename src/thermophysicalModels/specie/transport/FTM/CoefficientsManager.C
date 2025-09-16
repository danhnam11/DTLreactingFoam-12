/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::CoefficientsManager

Description


SourceFiles
    CoefficientsManager.H
    CoefficientsManager.C

\*---------------------------------------------------------------------------*/

#include "CoefficientsManager.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::CoefficientsManager::CoefficientsManager()
:
    muCoeffs_(),
    kappaCoeffs_(),
    DijCoeffs_()
{
}

Foam::CoefficientsManager::CoefficientsManager(const Foam::dictionary& dict)
:
    muCoeffs_(),
    kappaCoeffs_(),
    DijCoeffs_()
{
    dict.lookup("muCoeffs") >> muCoeffs_;
    dict.lookup("kappaCoeffs") >> kappaCoeffs_;
    dict.lookup("DijCoeffs") >> DijCoeffs_;

    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CoefficientsManager::setMkCoeffs
(
    const List<List<scalar>>& muCoeffsMk,
    const List<List<scalar>>& kappaCoeffsMk,
    const List<List<List<scalar>>>& DijCoeffsMk
)
{
    muCoeffsMk_ = muCoeffsMk;
    kappaCoeffsMk_ = kappaCoeffsMk;
    DijCoeffsMk_ = DijCoeffsMk;
}

void Foam::CoefficientsManager::setMuKappaMkCoeffs
(
    const List<List<scalar>>& muCoeffsMk,
    const List<List<scalar>>& kappaCoeffsMk
)
{
    muCoeffsMk_ = muCoeffsMk;
    kappaCoeffsMk_ = kappaCoeffsMk;
}


void Foam::CoefficientsManager::printCoeffs() const
{
   /*
    Info<< ">>> DEBUG: CoefficieintsManger constructed" << nl;
    Info << "muCoeffs_: " << muCoeffs_ << nl;
    Info << "kappaCoeffs_: " << kappaCoeffs_ << nl;
    Info << "DijCoeffs_: " << DijCoeffs_ << nl;

    Info << "muCoeffsMk_: " << muCoeffsMk_ << nl;
    Info << "kappaCoeffsMk_: " << kappaCoeffsMk_ << nl;
    Info << "DijCoeffsMk_: " << DijCoeffsMk_ << nl;
    */
}

const Foam::List<Foam::scalar>& Foam::CoefficientsManager::muCoeffs() const
{
    return muCoeffs_; 
}

const Foam::List<Foam::scalar>& Foam::CoefficientsManager::kappaCoeffs() const
{ 
    return kappaCoeffs_; 
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::CoefficientsManager::DijCoeffs() const
{ 
    return DijCoeffs_; 
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::CoefficientsManager::muCoeffsMk() const
{
    return muCoeffsMk_;
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::CoefficientsManager::kappaCoeffsMk() const
{
    return kappaCoeffsMk_;
}

const Foam::List<Foam::List<Foam::List<Foam::scalar>>>& Foam::CoefficientsManager::DijCoeffsMk() const
{
    return DijCoeffsMk_;
}


// ************************************************************************* //
