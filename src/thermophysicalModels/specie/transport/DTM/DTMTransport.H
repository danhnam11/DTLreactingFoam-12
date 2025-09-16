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
    Foam::DTMTransport

Description
    Transport package using Standard Kinetic Theory model. The mixture is 
    calculated based on the modification of the Wilke semi-empirical formulas. 

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
            // input parameters for DTM
            linearity       1;
            epsilonOverKb   107.4;
            sigma           3.458;
            dipoleMoment    0;
            alpha           1.6;
            Zrot            3.8;
        }
    }
    \endverbatim

SourceFiles
    DTMTransportI.H
    DTMTransport.C

\*---------------------------------------------------------------------------*/

#ifndef DTMTransport_H
#define DTMTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class DTMTransport;

template<class Thermo>
inline DTMTransport<Thermo> operator+
(
    const DTMTransport<Thermo>&,
    const DTMTransport<Thermo>&
);

template<class Thermo>
inline DTMTransport<Thermo> operator*
(
    const scalar,
    const DTMTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const DTMTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                     Class DTMTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class DTMTransport
:
    public Thermo
{
    // Private Data

protected:

     //- Members for diffusivity evaluation
        scalar linearity_;     //linearity
        scalar epsilonOverKb_; //Lennard-Jones well depth
        scalar sigma_;         //Lennard-Jones diameter
        scalar miui_;          //Dipole moment 
        scalar polar_;         //Polarizability
        scalar Zrot_;          //Rotational relaxation number
        scalar Kb_;

        //- list of mass fraction 
        List<scalar> Ymd_;

        //- list of mole fraction 
        List<scalar> Xmd_;

        //- list of epsilonijOverKb[i][j] 
        List<List<scalar>> epsilonijOverKb_;

        //- list of deltaij[i][j] 
        List<List<scalar>> deltaij_;

        //- list of Mij[i][j] 
        List<List<scalar>> Mij_;

        //- list of sigmaij[i][j] 
        List<List<scalar>> sigmaij_;

      //for visc. and cond.
        //- list of all members of k-th specie in the mixture
        List<scalar> linearityMk_;
        List<scalar> epsilonOverKbMk_;
        List<scalar> sigmaMk_;
        List<scalar> miuiMk_;
        List<scalar> polarMk_;
        List<scalar> ZrotMk_;
        List<scalar> wMk_;
        List<List<List<scalar>>> CpCoeffTableMk_;

        List<List<scalar>> muCoeffsMk_;
        List<List<scalar>> kappaCoeffsMk_;
        List<List<List<scalar>>> DijCoeffsMk_;

      //- protected functions
        //- Quadratic interpolation 
        inline scalar qInterpolate(scalar t, List<scalar> X, List<scalar> Y) const;

        // Return collision integral from tablulated data
        //- (Table IV, Monchick and Mason, Journal of Chem. Phys. 35: 1676 (1961))
        inline scalar OMEGA(int N, scalar TR, scalar DR) const;


        //- Funtion return Zrot(T)
        inline scalar ZrotT(scalar Zrot0, scalar epsilonOverKb, scalar T) const;

        //- Temporarily calculate the dynamic viscosity [kg/ms] of pure species
        inline scalar etak
        (
            const scalar Wk,
            const scalar epsilonOverKbk,
            const scalar sigmak,
            const scalar miuik,
            const scalar p, 
            const scalar T
        ) const;

        //- Temporarily calculate the thermal conductivity [W/mK] of pure species
        inline scalar lamdak
        (
            const scalar Wk,
            const scalar linearityk,
            const scalar epsilonOverKbk,
            const scalar sigmak,
            const scalar miuik,
            const scalar Zrotk,
            const List<List<scalar>> CpCoeffTablek,
            const scalar p, 
            const scalar T
        ) const;



public:

    // Constructors

        //- Construct from components
        inline DTMTransport
        (
            const Thermo& t,
            const scalar& linearity,
            const scalar& epsilonOverKb,
            const scalar& sigma,
            const scalar& miui,
            const scalar& polar,
            const scalar& Zrot,
            const scalar& Kb,

            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij,

            const List<scalar>& linearityMk,
            const List<scalar>& epsilonOverKbMk,
            const List<scalar>& sigmaMk,
            const List<scalar>& miuiMk,
            const List<scalar>& polarMk,
            const List<scalar>& ZrotMk,
            const List<scalar>& wMk,
            const List<List<List<scalar>>>& CpCoeffTableMk
        );


        //- Construct as named copy
        inline DTMTransport(const word&, const DTMTransport&);

        //- Construct from name and dictionary
        DTMTransport(const word& name, const dictionary& dict);

        //- Construct from base thermo and dictionary
        DTMTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<DTMTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "DTMTransport<" + Thermo::typeName() + '>';
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

        //- Update coefficients for mixture transport properties
        inline void updateTRANS
        (
            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij,
            
            const List<scalar>& linearityMk,
            const List<scalar>& epsilonOverKbMk,
            const List<scalar>& sigmaMk,
            const List<scalar>& miuiMk,
            const List<scalar>& polarMk,
            const List<scalar>& ZrotMk,
            const List<scalar>& wMk,
            const List<List<List<scalar>>>& CpCoeffTableMk 
        );

        //- Functions return member data for testing mixture  
        // for diffusivity
        inline scalar linearity() const;
        inline scalar epsilonOverKb() const;
        inline scalar sigma() const;
        inline scalar miui() const;
        inline scalar polar() const;
        inline scalar Zrot() const;
        inline scalar Kb() const;

        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;
        inline List<List<scalar>> epsilonijOverKb() const;
        inline List<List<scalar>> deltaij() const;
        inline List<List<scalar>> Mij() const;
        inline List<List<scalar>> sigmaij() const;

        inline List<List<scalar>> muCoeffsMk() const;
        inline List<List<scalar>> kappaCoeffsMk() const;
        inline List<List<List<scalar>>> DijCoeffsMk() const;


    // Member Operators

        inline void operator+=(const DTMTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend DTMTransport operator+ <Thermo>
        (
            const DTMTransport&,
            const DTMTransport&
        );

        friend DTMTransport operator* <Thermo>
        (
            const scalar,
            const DTMTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const DTMTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DTMTransportI.H"

#ifdef NoRepository
    #include "DTMTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
