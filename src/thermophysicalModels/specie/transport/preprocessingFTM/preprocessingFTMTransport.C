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

#include "preprocessingFTMTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- construct from components
template<class Thermo>
inline Foam::preprocessingFTMTransport<Thermo>::preprocessingFTMTransport
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
)
:
    Thermo(t),
    // for diffusivity
    linearity_(linearity),
    epsilonOverKb_(epsilonOverKb),
    sigma_(sigma),
    miui_(miui),
    polar_(polar),
    Zrot_(Zrot),
    Kb_(Kb),

    Ymd_(Ymd),
    Xmd_(Xmd),
    epsilonijOverKb_(epsilonijOverKb),
    deltaij_(deltaij),
    Mij_(Mij),
    sigmaij_(sigmaij),

    linearityMk_(linearityMk),
    epsilonOverKbMk_(epsilonOverKbMk),
    sigmaMk_(sigmaMk),
    miuiMk_(miuiMk),
    polarMk_(polarMk),
    ZrotMk_(ZrotMk),
    wMk_(wMk),
    CpCoeffTableMk_(CpCoeffTableMk)
{}


// construct as named copy
template<class Thermo>
inline Foam::preprocessingFTMTransport<Thermo>::preprocessingFTMTransport
(
    const word& name,
    const preprocessingFTMTransport& st
)
:
    Thermo(name, st),
    // for diffusivity
    linearity_(st.linearity_),
    epsilonOverKb_(st.epsilonOverKb_),
    sigma_(st.sigma_),
    miui_(st.miui_),
    polar_(st.polar_),
    Zrot_(st.Zrot_),
    Kb_(st.Kb_),    

    Ymd_(st.Ymd_),
    Xmd_(st.Xmd_),
    epsilonijOverKb_(st.epsilonijOverKb_),
    deltaij_(st.deltaij_),
    Mij_(st.Mij_),
    sigmaij_(st.sigmaij_),

    linearityMk_(st.linearityMk_),
    epsilonOverKbMk_(st.epsilonOverKbMk_),
    sigmaMk_(st.sigmaMk_),
    miuiMk_(st.miuiMk_),
    polarMk_(st.polarMk_),
    ZrotMk_(st.ZrotMk_),
    wMk_(st.wMk_),
    CpCoeffTableMk_(st.CpCoeffTableMk_)
{}


//- Construct from name and dictionary
template<class Thermo>
Foam::preprocessingFTMTransport<Thermo>::preprocessingFTMTransport
(
    const word& name,
    const dictionary& dict
)
:
    Thermo(name, dict),
    // for diffusivity
    linearity_(dict.subDict("transport").lookup<scalar>("linearity")),
    epsilonOverKb_(dict.subDict("transport").lookup<scalar>("epsilonOverKb")),
    sigma_(dict.subDict("transport").lookup<scalar>("sigma")),
    miui_(dict.subDict("transport").lookup<scalar>("dipoleMoment")),
    polar_(dict.subDict("transport").lookup<scalar>("alpha")),
    Zrot_(dict.subDict("transport").lookup<scalar>("Zrot")),  
    Kb_(8.314510/(6.0221367*1E23)),
 
    Ymd_(2),
    Xmd_(2),
    epsilonijOverKb_(2),
    deltaij_(2),
    Mij_(2),
    sigmaij_(2),

    linearityMk_(2),
    epsilonOverKbMk_(2),
    sigmaMk_(2),
    miuiMk_(2),
    polarMk_(2),
    ZrotMk_(2),
    wMk_(2),
    CpCoeffTableMk_(2)
{
  //- Temporary initialization
    forAll(Ymd_, i) 
    {
        //Ymd_[i] = this->Y();
        //Xmd_[i] = this->Y()/this->W();

        Ymd_[i] = 0.5;
        Xmd_[i] = 0.5;

        linearityMk_[i] = linearity_;
        epsilonOverKbMk_[i] = epsilonOverKb_;
        sigmaMk_[i] = sigma_;
        miuiMk_[i] = miui_;
        polarMk_[i] = polar_;
        ZrotMk_[i] = Zrot_;
        wMk_[i] = this->W();
        CpCoeffTableMk_[i] = this->CpCoeffTable();
    } 

    forAll(epsilonijOverKb_, i) 
    {
        epsilonijOverKb_[i] = epsilonOverKb_;
        deltaij_[i] = 1.0;
        Mij_[i] = this->W();
        sigmaij_[i] = sigma_;
    }
}



//- Construct from base thermo and dictionary
template<class Thermo>
Foam::preprocessingFTMTransport<Thermo>::preprocessingFTMTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    // for diffusivity
    linearity_(dict.subDict("transport").lookup<scalar>("linearity")),
    epsilonOverKb_(dict.subDict("transport").lookup<scalar>("epsilonOverKb")),
    sigma_(dict.subDict("transport").lookup<scalar>("sigma")),
    miui_(dict.subDict("transport").lookup<scalar>("dipoleMoment")),
    polar_(dict.subDict("transport").lookup<scalar>("alpha")),
    Zrot_(dict.subDict("transport").lookup<scalar>("Zrot")),  
    Kb_(8.314510/(6.0221367*1E23)),
 
    Ymd_(2),
    Xmd_(2),
    epsilonijOverKb_(2),
    deltaij_(2),
    Mij_(2),
    sigmaij_(2),

    linearityMk_(2),
    epsilonOverKbMk_(2),
    sigmaMk_(2),
    miuiMk_(2),
    polarMk_(2),
    ZrotMk_(2),
    wMk_(2),
    CpCoeffTableMk_(2)
{
  //- Temporary initialization
    forAll(Ymd_, i) 
    {
        Ymd_[i] = this->Y();
        Xmd_[i] = this->Y()/this->W();

        linearityMk_[i] = linearity_;
        epsilonOverKbMk_[i] = epsilonOverKb_;
        sigmaMk_[i] = sigma_;
        miuiMk_[i] = miui_;
        polarMk_[i] = polar_;
        ZrotMk_[i] = Zrot_;
        wMk_[i] = this->W();
        CpCoeffTableMk_[i] = this->CpCoeffTable();
    } 

    forAll(epsilonijOverKb_, i) 
    {
        epsilonijOverKb_[i] = epsilonOverKb_;
        deltaij_[i] = 1.0;
        Mij_[i] = this->W();
        sigmaij_[i] = sigma_;
    }
}



//- Construct and return a clone
template<class Thermo>
inline Foam::autoPtr<Foam::preprocessingFTMTransport<Thermo>>
Foam::preprocessingFTMTransport<Thermo>::clone() const
{
    return autoPtr<preprocessingFTMTransport<Thermo>>
    (
        new preprocessingFTMTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::preprocessingFTMTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("linearity", linearity_);
    dict.add("epsilonOverKb", epsilonOverKb_);
    dict.add("sigma", sigma_);
    dict.add("dipoleMoment", miui_);
    dict.add("alpha", polar_);
    dict.add("Zrot", Zrot_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const preprocessingFTMTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::preprocessingFTMTransport<Thermo>::operator+=
(
    const preprocessingFTMTransport<Thermo>& st
)
{

    Thermo::operator+=(st);

        // for diffusivity
        linearity_       = st.linearity_;
        epsilonOverKb_   = st.epsilonOverKb_;
        sigma_           = st.sigma_;
        miui_            = st.miui_;
        polar_           = st.polar_;
        Zrot_            = st.Zrot_;
        Kb_              = st.Kb_;

        Ymd_             = st.Ymd_; 
        Xmd_             = st.Xmd_; 
        epsilonijOverKb_ = st.epsilonijOverKb_; 
        deltaij_         = st.deltaij_; 
        Mij_             = st.Mij_; 
        sigmaij_         = st.sigmaij_;

        linearityMk_       = st.linearityMk_;
        epsilonOverKbMk_   = st.epsilonOverKbMk_;
        sigmaMk_           = st.sigmaMk_;
        miuiMk_            = st.miuiMk_;
        polarMk_           = st.polarMk_;
        ZrotMk_            = st.ZrotMk_;
        wMk_               = st.wMk_;
        CpCoeffTableMk_    = st.CpCoeffTableMk_;
}


template<class Thermo>
inline void Foam::preprocessingFTMTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::preprocessingFTMTransport<Thermo> Foam::operator+
(
    const preprocessingFTMTransport<Thermo>& st1,
    const preprocessingFTMTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < small)
    {
        return preprocessingFTMTransport<Thermo>
        (
            t,
            // for diffusivity
            st1.linearity_,
            st1.epsilonOverKb_,
            st1.sigma_,
            st1.miui_,
            st1.polar_,
            st1.Zrot_,
            st1.Kb_,

            st1.Ymd_,
            st1.Xmd_,
            st1.epsilonijOverKb_,
            st1.deltaij_,
            st1.Mij_,
            st1.sigmaij_,

            st1.linearityMk_,
            st1.epsilonOverKbMk_,
            st1.sigmaMk_,
            st1.miuiMk_,
            st1.polarMk_,
            st1.ZrotMk_,
            st1.wMk_,
            st1.CpCoeffTableMk_
        );
    }
    else
    {
        //scalar Y1 = st1.Y()/t.Y();
        //scalar Y2 = st2.Y()/t.Y();

        return preprocessingFTMTransport<Thermo>
        (
            t,
            // for diffusivity
            st1.linearity_,
            st1.epsilonOverKb_,
            st1.sigma_,
            st1.miui_,
            st1.polar_,
            st1.Zrot_,
            st1.Kb_,

            st1.Ymd_,
            st1.Xmd_,
            st1.epsilonijOverKb_,
            st1.deltaij_,
            st1.Mij_,
            st1.sigmaij_,

            st1.linearityMk_,
            st1.epsilonOverKbMk_,
            st1.sigmaMk_,
            st1.miuiMk_,
            st1.polarMk_,
            st1.ZrotMk_,
            st1.wMk_,
            st1.CpCoeffTableMk_
        );
    }
}


template<class Thermo>
inline Foam::preprocessingFTMTransport<Thermo> Foam::operator*
(
    const scalar s,
    const preprocessingFTMTransport<Thermo>& st
)
{
    return preprocessingFTMTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        // for diffusivity
        st.linearity_,
        st.epsilonOverKb_,
        st.sigma_,
        st.miui_,
        st.polar_,
        st.Zrot_,
        st.Kb_,

        st.Ymd_,
        st.Xmd_,
        st.epsilonijOverKb_,
        st.deltaij_,
        st.Mij_,
        st.sigmaij_,

        st.linearityMk_,
        st.epsilonOverKbMk_,
        st.sigmaMk_,
        st.miuiMk_,
        st.polarMk_,
        st.ZrotMk_,
        st.wMk_,
        st.CpCoeffTableMk_
    );
}


// ************************************************************************* //
