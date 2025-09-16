/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#ifndef PROPERTYREADER_H
#define PROPERTYREADER_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "fvm.H"

#include "List.H"
#include "scalar.H"

//#include <sys/stat.h>
//#include <sys/types.h>
//#include <iostream>

#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace Foam;

class PropertyReader
{
public:
    PropertyReader
    (
        const Foam::fvMesh& mesh,
	const Foam::wordList& speciesList
    );

    void generateTranCoeff(const Foam::wordList& speciesList);
    void printTranCoeff
    (
        const Foam::Time& runTime, 
	const Foam::word& species,
	List<scalar>& muTci,
	List<scalar>& kappaTci,
	List<List<scalar>>& DijTci
    );


private:
    Foam::wordList propertiesList;
    Foam::PtrList<Foam::volScalarField> fields;
    Foam::List<Foam::List<Foam::volScalarField*>> dijFields;
    Foam::word caseDirectory;
    const Foam::fvMesh& mesh;
    Foam::wordList speciesNames;
};

#endif // PROPERTYREADER_H
