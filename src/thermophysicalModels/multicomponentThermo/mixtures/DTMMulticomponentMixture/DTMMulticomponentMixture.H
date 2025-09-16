/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::DTMMulticomponentMixture

Description
    Thermophysical properties mixing class which applies mass-fraction weighted
    mixing to the thermodynamic coefficients.

    Transport properties are based on Standard Kinetic Theory model. 
    The mixture is calculated based on the modification of the Wilke semi-empirical formulas. 

    Mass Diffusivity is based on mixture averaged model in which binary diffusion 
    coefficients are obtained from Standard Kinetic Theory model.

    It has been validated against NIST, CHEMKIN, Cantera, laminarSMOKE
 
    For details of equations in the model, refer to TRANSPORT program 
    in CHEMKIN II @ R. J. Kee et al. CHEMKIN collection, Release 3.6, 
    Reaction Design, Inc., San Diego, CA (2001). 

    Output: 
        Dynamic Viscosity [kg/m.s]
        Thermal Conductivity [W/m.K]    
        Mass Diffusivity [m2/s]    

    developed by:
    Nam Danh Nguyen and Jae Hun Lee
    Advisor: Prof. Chun Sang Yoo 
    Clean combustion & energy research (CCER) Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    DTMMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef DTMMulticomponentMixture_H
#define DTMMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class DTMMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class DTMMulticomponentMixture
:
    public multicomponentMixture<ThermoType>
{
public:

    // Public Typedefs

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;

	//- for FTM
        typedef DTMMulticomponentMixture<ThermoType> mixtureType;
        typedef ThermoType thermoType;

private:

    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable thermoMixtureType mixture_;


        //- For mixture DTM - Nam
        //- List of species molecular weight
        List<scalar> ListW_;

        List<scalar> linearityMk_;
        List<scalar> epsilonOverKbMk_;
        List<scalar> sigmaMk_;
        List<scalar> miuiMk_;
        List<scalar> polarMk_;
        List<scalar> ZrotMk_;
        List<List<List<scalar>>> CpCoeffTableMk_;

        List<List<scalar>> EPSILONijOVERKB_;
        List<List<scalar>> DELTAij_;
        List<List<scalar>> Mij_;
        List<List<scalar>> SIGMAij_;
        // 



public:

    // Constructors

        //- Construct from a dictionary
        DTMMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        DTMMulticomponentMixture
        (
            const DTMMulticomponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~DTMMulticomponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "DTMMulticomponentMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;

        inline const PtrList<ThermoType>& specieThermos() const
        {
            return multicomponentMixture<ThermoType>::specieThermos();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DTMMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
