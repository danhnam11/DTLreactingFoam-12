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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// moved to .C file
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar Eta = 0;

    forAll(Xmd_, k)
    {
        scalar Etak = etak(k, wMk_[k], p, T); 
        scalar sumPhi = 0;

        forAll(Xmd_, j)
        {
            scalar WkOj = wMk_[k]/wMk_[j];
            scalar Etaj = etak(j, wMk_[j], p, T);
            scalar X = 1 + sqrt(Etak/Etaj) * pow(wMk_[j]/wMk_[k], 0.25);

            sumPhi = sumPhi + Xmd_[j]*(pow((1 + WkOj), -0.5) * pow(X, 2.0)/sqrt(8.0));
        }

        //Eta = Eta + Xmd_[k]*Etak/(sumPhi+1e-24);
        Eta = Eta + Xmd_[k]*Etak/(sumPhi);

    }
    if (Eta < 0) {Eta =  1e-16;}  
    //return (0.1*Eta); // for DTM
    return (Eta); // for FTM
}


template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    scalar sum1 = 0; 
    scalar sum2 = 0;
    forAll(Xmd_, k)
    {
        scalar Lamdak = lamdak(k, wMk_[k], p, T);
        sum1 = sum1 + Xmd_[k]*Lamdak;
        sum2 = sum2 + Xmd_[k]/Lamdak;
    }

    scalar Lamda = 0.5*(sum1 + 1/(sum2 + 1e-16));
    return Lamda; //in si unit
}


// mass Diffusivity 
template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::Dimix
(
    label speciei,
    const scalar p,
    const scalar T
) const
{
// for FTM

    scalar sumXW  = 0.0;
    scalar sumXOD = 0.0;
    scalar WTM    = 0.0;
    scalar Dimix  = 0.0;

    forAll(wMk_,j)
    {
        WTM = WTM + wMk_[j]*Xmd_[j];
    }

    forAll(Xmd_, j)
    {
        if (speciei != j)
        {
            sumXW  = sumXW + Xmd_[j] * wMk_[j];
            sumXOD += Xmd_[j] / (DijFit(speciei, j, p, T));
        }
    }
    Dimix = sumXW/(WTM*sumXOD);

    return Dimix;
} 


// thermal diffusivity
template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::DimixT
(
    label speciei,
    const scalar p,
    const scalar T
) const
{
         scalar sumT1 = 0.0;
         scalar sumT2 = 0.0;
         scalar DimixT = 0.0;

         forAll(Xmd_, j)
         {
             sumT1 = sumT1 + pow(wMk_[j],0.511)*Xmd_[j];
             sumT2 = sumT2 + pow(wMk_[j],0.489)*Xmd_[j];
         }
         DimixT = (-2.59*1e-7*pow(T,0.659))*(pow(wMk_[speciei],0.511)*Xmd_[speciei]/sumT1-Ymd_[speciei])*(sumT1/sumT2);

     return DimixT;
}

template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::DijFit
(
    label speciei,
    label speciej,
    const scalar p,
    const scalar T
) const
{

    scalar DijFit = 0;

    if (speciei == speciej)
    {
        DijFit = 1.0;
    }
    else
    {
        scalar aD0 = this->coeffs_.DijCoeffsMk()[speciei][speciej][0];
        scalar aD1 = this->coeffs_.DijCoeffsMk()[speciei][speciej][1];
        scalar aD2 = this->coeffs_.DijCoeffsMk()[speciei][speciej][2];
        scalar aD3 = this->coeffs_.DijCoeffsMk()[speciei][speciej][3];

        scalar lnT = log(T);
        DijFit = exp(aD0 + aD1*lnT+ aD2*pow(lnT,2) + aD3 * pow(lnT,3));
    }
    return DijFit; // cgs to SI unit

}

// update coefficients for mixture calculation
template<class Thermo>
inline void Foam::FTMTransport<Thermo>::updateTRANS
(
    const List<scalar>& Ymd,
    const List<scalar>& Xmd,
    const List<scalar>& wMk
)
{
    Ymd_             = Ymd;
    Xmd_             = Xmd;
    wMk_             = wMk;

}

template<class Thermo>
inline void Foam::FTMTransport<Thermo>::updateTRANSFitCoeff
(
    const List<List<scalar>>& muCoeffsMk,
    const List<List<scalar>>& kappaCoeffsMk,
    const List<List<List<scalar>>>& DijCoeffsMk
)
{
    coeffs_.setMkCoeffs(muCoeffsMk, kappaCoeffsMk, DijCoeffsMk);
}

// functions to access coefficients
template<class Thermo>
inline Foam::List<Foam::scalar>
Foam::FTMTransport<Thermo>::Ymd() const
{
    return Ymd_;
}

template<class Thermo>
inline Foam::List<Foam::scalar>
Foam::FTMTransport<Thermo>::Xmd() const
{
    return Xmd_;
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

//- Temporarily calculate the dynamic viscosity of pure species
template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::etak
(
    label speciei,
    const scalar Wk,
    const scalar p,
    const scalar T
) const
{
    scalar etak = 0.0;

    scalar am0 = this->coeffs_.muCoeffsMk()[speciei][0];
    scalar am1 = this->coeffs_.muCoeffsMk()[speciei][1];
    scalar am2 = this->coeffs_.muCoeffsMk()[speciei][2];
    scalar am3 = this->coeffs_.muCoeffsMk()[speciei][3];
    scalar lnT = log(T);
    etak = exp(am0 + am1*lnT+ am2*pow(lnT,2) + am3 * pow(lnT,3));

    return etak; // in cgs unist

}


//- Temporarily calculate the thermal conductivity [W/mK] of pure species
template<class Thermo>
inline Foam::scalar Foam::FTMTransport<Thermo>::lamdak
(
    label speciei,
    const scalar Wk,
    const scalar p,
    const scalar T
) const
{
    scalar lamdak = 0.0;

    scalar ak0 = this->coeffs_.kappaCoeffsMk()[speciei][0];
    scalar ak1 = this->coeffs_.kappaCoeffsMk()[speciei][1];
    scalar ak2 = this->coeffs_.kappaCoeffsMk()[speciei][2];
    scalar ak3 = this->coeffs_.kappaCoeffsMk()[speciei][3];

    scalar lnT = log(T);
    lamdak = exp(ak0 + ak1*lnT+ ak2*pow(lnT,2) + ak3 * pow(lnT,3));

    return lamdak;
}


// ************************************************************************* //
