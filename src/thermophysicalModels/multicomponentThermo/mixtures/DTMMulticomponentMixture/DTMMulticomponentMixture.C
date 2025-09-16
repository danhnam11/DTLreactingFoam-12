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

#include "DTMMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::DTMMulticomponentMixture<ThermoType>::
DTMMulticomponentMixture
(
    const dictionary& dict
)
:
    multicomponentMixture<ThermoType>(dict),
    mixture_("mixture", this->specieThermos()[0]),

    // for mixture DTM - Nam
    ListW_(this->specieThermos().size()),
    linearityMk_(this->specieThermos().size()),
    epsilonOverKbMk_(this->specieThermos().size()),
    sigmaMk_(this->specieThermos().size()),
    miuiMk_(this->specieThermos().size()),
    polarMk_(this->specieThermos().size()),
    ZrotMk_(this->specieThermos().size()),
    CpCoeffTableMk_(this->specieThermos().size()),

    EPSILONijOVERKB_(this->specieThermos().size()),
    DELTAij_(this->specieThermos().size()),
    Mij_(this->specieThermos().size()),
    SIGMAij_(this->specieThermos().size())
    //
{

    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]           = this->specieThermos()[i].W();
        linearityMk_[i]     = this->specieThermos()[i].linearity();
        epsilonOverKbMk_[i] = this->specieThermos()[i].epsilonOverKb();
        sigmaMk_[i]         = this->specieThermos()[i].sigma();
        miuiMk_[i]          = this->specieThermos()[i].miui();
        polarMk_[i]         = this->specieThermos()[i].polar();
        ZrotMk_[i]          = this->specieThermos()[i].Zrot();        
        CpCoeffTableMk_[i]  = this->specieThermos()[i].CpCoeffTable();
    }

    List<scalar> nEPSILONijOVERKB(this->specieThermos().size());
    List<scalar> nDELTAij(this->specieThermos().size());
    List<scalar> nMij(this->specieThermos().size());
    List<scalar> nSIGMAij(this->specieThermos().size());    

    // Dip = dipole moment 
    const scalar DipMin = 1e-20;

    forAll(EPSILONijOVERKB_, i)
    {
        forAll(nEPSILONijOVERKB, j)
        {
            nMij[j] = 
               1/(1/this->specieThermos()[i].W() + 1/this->specieThermos()[j].W());

            nEPSILONijOVERKB[j] = 
               sqrt(this->specieThermos()[i].epsilonOverKb()*this->specieThermos()[j].epsilonOverKb());

            nSIGMAij[j] = 
               0.5*(this->specieThermos()[i].sigma() + this->specieThermos()[j].sigma());

            // calculate coeficient xi
            scalar xi = 1.0; 
            if ((this->specieThermos()[i].miui() < DipMin) && (this->specieThermos()[j].miui() > DipMin)) 
            {
                // miui_j > DipMin > miui_i --> j is polar, i is nonpolar              
                nDELTAij[j] = 0;
                xi = 1.0 + 
                     this->specieThermos()[i].polar()*pow(this->specieThermos()[j].miui(), 2)* 
                     sqrt(this->specieThermos()[j].epsilonOverKb()/this->specieThermos()[i].epsilonOverKb())/
                     (4*pow(this->specieThermos()[i].sigma(), 3)*
                            (
                             1e+19*this->specieThermos()[j].Kb()*this->specieThermos()[j].epsilonOverKb()*
                             pow(this->specieThermos()[j].sigma(), 3)
                            )
                     );
            }
            else if ((this->specieThermos()[i].miui() > DipMin) && (this->specieThermos()[j].miui() < DipMin))
            {
                // miui_j < DipMin < miui_i --> i is polar, j is nonpolar
                nDELTAij[j] = 0;
                xi = 1.0 +
                     this->specieThermos()[j].polar()*pow(this->specieThermos()[i].miui(), 2)*
                     sqrt(this->specieThermos()[i].epsilonOverKb()/this->specieThermos()[j].epsilonOverKb())/
                     (4*pow(this->specieThermos()[j].sigma(), 3)*   
                            (
                             1e+19*this->specieThermos()[i].Kb()*this->specieThermos()[i].epsilonOverKb()*
                             pow(this->specieThermos()[i].sigma(), 3)
                            )
                     ); 
            }
            else 
            {            
                xi = 1.0;
                nDELTAij[j] =
                   0.5*(this->specieThermos()[i].miui()*this->specieThermos()[j].miui())/
                   (nEPSILONijOVERKB[j]*1e+19*this->specieThermos()[j].Kb()*pow(nSIGMAij[j], 3));
            }
           
            nEPSILONijOVERKB[j] = nEPSILONijOVERKB[j]*pow(xi, 2);
            nSIGMAij[j] = nSIGMAij[j]*pow(xi, -1/6);
        }

        EPSILONijOVERKB_[i] = nEPSILONijOVERKB;
        Mij_[i]             = nMij;
        SIGMAij_[i]         = nSIGMAij;
        DELTAij_[i]         = nDELTAij;
    }
    // End of pre-calculation for kinetic model

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::DTMMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::DTMMulticomponentMixture<ThermoType>::thermoMixture
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
    mixture_.updateTRANS
    (
        listYi, listXi, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
    // 

    return mixture_;
}


template<class ThermoType>
const typename
Foam::DTMMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::DTMMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return thermoMixture(Y);
}


template<class ThermoType>
const typename
Foam::DTMMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::DTMMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}

// ************************************************************************* //
