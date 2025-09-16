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
    Foam::FTMTransport

Description
    Transport package using Fitting Polynomial Model. 
    Coefficients are generated based on Standard Kinetic Theory model. 
    The mixture is calculated based on the modification of the Wilke 
    semi-empirical formulas. 

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

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

Usage

    Example specification of DTMTransport for gases:
    \verbatim
    O2
    {
        ...
        transport
        {
            // input parameters for FTM
            muCoeffs        (-19.3405 2.62765 -0.265133 0.0118167);
            kappaCoeffs     (-13.8944 3.07851 -0.297333 0.0127535);
            DijCoeffs       (
                                (-24.3747 3.3315 -0.221273 0.00973849)
                                (-30.453 5.60295 -0.492151 0.0205307)
                                (-9.21034 -6.45e-26 9.26e-27 -4.40e-28)
                                (-27.6164 4.47195 -0.36524 0.0158195)
                                (-25.3927 3.75885 -0.275697 0.0120551)
                                ...
                            );
        }
    }
    \endverbatim

SourceFiles
    FTMTransportI.H
    FTMTransport.C

\*---------------------------------------------------------------------------*/

#ifndef FTMTransport_H
#define FTMTransport_H

#include "CoefficientsManager.H" //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class FTMTransport;

template<class Thermo>
inline FTMTransport<Thermo> operator+
(
    const FTMTransport<Thermo>&,
    const FTMTransport<Thermo>&
);

template<class Thermo>
inline FTMTransport<Thermo> operator*
(
    const scalar,
    const FTMTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const FTMTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                     Class FTMTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class FTMTransport
:
    public Thermo
{
    // Private Data
    CoefficientsManager internalCoeffs_;
    CoefficientsManager& coeffs_;

protected:

     //- Members for diffusivity evaluation
        //- list of mass fraction 
        List<scalar> Ymd_;

        //- list of mole fraction 
        List<scalar> Xmd_;

        //- list of all members of k-th specie in the mixture
        List<scalar> wMk_;

      //- protected functions
        //- Temporarily calculate the dynamic viscosity [kg/ms] of pure species
        inline scalar etak
        (
            label speciei,
            const scalar Wk,
            const scalar p, 
            const scalar T
        ) const;

        //- Temporarily calculate the thermal conductivity [W/mK] of pure species
        inline scalar lamdak
        (
            label speciei,
            const scalar Wk,
            const scalar p, 
            const scalar T
        ) const;

public:

    // Constructors

        //- Construct from components
        inline FTMTransport
        (
            const Thermo& t,
            CoefficientsManager& coeffs,

            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<scalar>& wMk
        );


        //- Construct as named copy
        inline FTMTransport(const word&, const FTMTransport&);

        //- Construct from name and dictionary
        FTMTransport(const word& name, const dictionary& dict);

        //- Construct from base thermo and dictionary
        FTMTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<FTMTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "FTMTransport<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;

        // 
        // Species diffusivity using Kinetic model for binary diffusion coefficients
        inline scalar Dimix(label speciei, const scalar p, const scalar T) const;

        // Species thermal diffusivity using Kinetic model for binary diffusion coefficients
        inline scalar DimixT(label speciei, const scalar p, const scalar T) const;

        // Binary diffusion coefficient
        inline scalar Dij(label speciei, label speciej, const scalar p, const scalar T) const
        {
            return 1.0;
        }

        // Binary diffusion coefficients Fit function
        inline scalar DijFit
        (
            label speciei,
            label speciej,
            const scalar p,
            const scalar T
        ) const;

        //- Update coefficients for mixture transport properties
        inline void updateTRANS
        (
            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<scalar>& wMk
        );

        //- Update coefficient for Fitting transport properties
        inline void updateTRANSFitCoeff
        (
            const List<List<scalar>>& muCoeffsMk,
            const List<List<scalar>>& kappaCoeffsMk,
            const List<List<List<scalar>>>& DijCoeffsMk
        );

        //- Functions return member data for testing mixture  
        // for diffusivity

        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;

        const CoefficientsManager& coeffs() const {return coeffs_;}

    // Member Operators

        inline void operator=(const FTMTransport&); //
        inline void operator+=(const FTMTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend FTMTransport operator+ <Thermo>
        (
            const FTMTransport&,
            const FTMTransport&
        );

        friend FTMTransport operator* <Thermo>
        (
            const scalar,
            const FTMTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const FTMTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FTMTransportI.H"

#ifdef NoRepository
    #include "FTMTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
