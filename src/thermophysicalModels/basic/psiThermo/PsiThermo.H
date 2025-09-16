/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::PsiThermo

Description
    Thermo implementation based on compressibility

SourceFiles
    PsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef PsiThermo_H
#define PsiThermo_H

#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class PsiThermo
:
    public BaseThermo
{

    // Private data 

        //- Store List of raw thermophysical data of all individual species - Nam      
        const PtrList<typename BaseThermo::mixtureType::thermoType>& specieThermos_;

    // Private Member Functions
        //- Calculate the thermo variables, original
        void calculate();

        // Nam 
        //- Initialize the thermo variables with full DTM/FTM
        void initialize();
        //- Calculate the thermo variables with coTHERM
        void calculateUsingCoTHERM();
        //- Calculate the thermo variables with coTHERM only T
        void calculateUsingCoTHERMOnlyT();
        //- Calculate the thermo variables, only for preprocessing FTM
        void calculateTransportPreProcessing();
        // Nam


public:

    // Constructors

        //- Construct from mesh and phase name
        PsiThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        PsiThermo(const PsiThermo<BaseThermo>&) = delete;


    //- Destructor
    virtual ~PsiThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PsiThermo<BaseThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
