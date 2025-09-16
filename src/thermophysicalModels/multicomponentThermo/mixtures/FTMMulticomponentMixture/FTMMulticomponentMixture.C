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

\*---------------------------------------------------------------------------*/

#include "FTMMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FTMMulticomponentMixture<ThermoType>::
FTMMulticomponentMixture
(
    const dictionary& dict
)
:
    multicomponentMixture<ThermoType>(dict),
    mixture_("mixture", this->specieThermos()[0]),

    // for mixture DTM - Nam
    ListW_(this->specieThermos().size()),
    muCoeffsMk_(this->specieThermos().size()),
    kappaCoeffsMk_(this->specieThermos().size()),
    DijCoeffsMk_(this->specieThermos().size())
    //
{

    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]         = this->specieThermos()[i].W();
        muCoeffsMk_[i]    = this->specieThermos()[i].coeffs().muCoeffs();
        kappaCoeffsMk_[i] = this->specieThermos()[i].coeffs().kappaCoeffs();
        DijCoeffsMk_[i]   = this->specieThermos()[i].coeffs().DijCoeffs();
    }
    // End of pre-calculation for kinetic model

    mixture_.updateTRANSFitCoeff
    (
        muCoeffsMk_, kappaCoeffsMk_, DijCoeffsMk_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::FTMMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::FTMMulticomponentMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    mixture_ = Y[0]*this->specieThermos()[0];

    for (label i=1; i<Y.size(); i++)
    {
        mixture_ += Y[i]*this->specieThermos()[i];
    }

    //- Update coefficients for mixture DTM - Nam
    //- List of secies mole and mass fraction 
    List<scalar> listXi(Y.size()); 
    List<scalar> listYi(Y.size()); 
    scalar sumXb = 0.0;  
    forAll(listXi, i)
    {
        sumXb = sumXb + Y[i]/ListW_[i]; 
    }  
    if (sumXb == 0){ sumXb = 1e-30;} 

    forAll(listXi, i)
    {
        listXi[i] = (Y[i]/ListW_[i])/sumXb;
        listYi[i] = Y[i];
        if(listXi[i] <= 0) { listXi[i] = 0; }
        if(listYi[i] <= 0) { listYi[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(listXi, i)
    {
        listXi[i] = listXi[i] + 1e-40;
        sumXcorrected = sumXcorrected + listXi[i];
    }
    
    forAll(listXi, i)
    {
        listXi[i] = listXi[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + listXi[i]*ListW_[i];
    }

    forAll(listYi, i)
    {
        listYi[i] = listXi[i]*ListW_[i]/WmixCorrect;
    }

    // Update coefficients for mixture DTM - Nam
    mixture_.updateTRANS(listYi, listXi, ListW_);
    //

    return mixture_;
}


template<class ThermoType>
const typename
Foam::FTMMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::FTMMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return thermoMixture(Y);
}


template<class ThermoType>
const typename
Foam::FTMMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::FTMMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


// ************************************************************************* //
