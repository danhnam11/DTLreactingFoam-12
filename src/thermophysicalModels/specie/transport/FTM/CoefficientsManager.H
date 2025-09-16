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

#ifndef CoefficientsManager_H
#define CoefficientsManager_H

#include "List.H"
#include "scalar.H"
#include "dictionary.H"
#include "IOstreams.H"

/*---------------------------------------------------------------------------*\
                 Class CoefficientsManager Declaration
\*---------------------------------------------------------------------------*/
namespace Foam
{

class CoefficientsManager
{
private:

    // Diffusivity evaluation coefficients
    List<scalar> muCoeffs_;
    List<scalar> kappaCoeffs_;
    List<List<scalar>> DijCoeffs_;

    List<List<scalar>> muCoeffsMk_;
    List<List<scalar>> kappaCoeffsMk_;
    List<List<List<scalar>>> DijCoeffsMk_;


public:
    // Default constructor
    CoefficientsManager();

    // Construct directly from dictionary
    explicit CoefficientsManager(const dictionary& dict);

    // Accessors
    const List<scalar>& muCoeffs() const;
    const List<scalar>& kappaCoeffs() const;
    const List<List<scalar>>& DijCoeffs() const;
    const List<List<scalar>>& muCoeffsMk() const;
    const List<List<scalar>>& kappaCoeffsMk() const;
    const List<List<List<scalar>>>& DijCoeffsMk() const;

    // Setter for Mk coefficients
    void setMkCoeffs
    (
        const List<List<scalar>>& muCoeffsMk,
        const List<List<scalar>>& kappaCoeffsMk,
        const List<List<List<scalar>>>& DijCoeffsMk
    );

    void setMuKappaMkCoeffs
    (
        const List<List<scalar>>& muCoeffsMk,
        const List<List<scalar>>& kappaCoeffsMk
    );

    // Optional: print function for debugging
    void printCoeffs() const;
    //void addSpecies(const dictionary& dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


