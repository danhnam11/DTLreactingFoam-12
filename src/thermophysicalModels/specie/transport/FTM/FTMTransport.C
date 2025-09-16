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

#include "FTMTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- construct from components
template<class Thermo>
inline Foam::FTMTransport<Thermo>::FTMTransport
(
    const Thermo& t,
    CoefficientsManager& coeffs,

    const List<scalar>& Ymd,
    const List<scalar>& Xmd,
    const List<scalar>& wMk
)
:
    Thermo(t),
    coeffs_(coeffs),

    // for diffusivity
    Ymd_(Ymd),
    Xmd_(Xmd),
    wMk_(wMk)
{}


// construct as named copy
template<class Thermo>
inline Foam::FTMTransport<Thermo>::FTMTransport
(
    const word& name,
    const FTMTransport& st
)
:
    Thermo(name, st),
    coeffs_(st.coeffs_),

    // for diffusivity
    Ymd_(st.Ymd_),
    Xmd_(st.Xmd_),
    wMk_(st.wMk_)
{
}


//- Construct from name and dictionary
template<class Thermo>
Foam::FTMTransport<Thermo>::FTMTransport
(
    const word& name,
    const dictionary& dict
)
:
    Thermo(name, dict),
    internalCoeffs_(dict.subDict("transport")),
    coeffs_(internalCoeffs_),
 
    // for diffusivity
 
    Ymd_(2),
    Xmd_(2),
    wMk_(2)
{
  //- Temporary initialization
    forAll(Ymd_, i) 
    {
        Ymd_[i] = this->Y();
        Xmd_[i] = this->Y()/this->W();
        wMk_[i] = this->W();
    } 
}



//- Construct from base thermo and dictionary
template<class Thermo>
Foam::FTMTransport<Thermo>::FTMTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    internalCoeffs_(dict.subDict("transport")),
    coeffs_(internalCoeffs_),

    // for diffusivity
    Ymd_(2),
    Xmd_(2),
    wMk_(2)
{
  //- Temporary initialization
    forAll(Ymd_, i) 
    {
        Ymd_[i] = this->Y();
        Xmd_[i] = this->Y()/this->W();
        wMk_[i] = this->W();
    } 
}

//- Construct and return a clone
template<class Thermo>
inline Foam::autoPtr<Foam::FTMTransport<Thermo>>
Foam::FTMTransport<Thermo>::clone() const
{
    return autoPtr<FTMTransport<Thermo>>
    (
        new FTMTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::FTMTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    //dict.add("As", As_);
    //dict.add("Ts", Ts_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const FTMTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::FTMTransport<Thermo>::operator=
(
    const FTMTransport<Thermo>& st
)
{
    Thermo::operator=(st); // JH modified

    Ymd_             = st.Ymd_;
    Xmd_             = st.Xmd_;
    wMk_             = st.wMk_;
}


template<class Thermo>
inline void Foam::FTMTransport<Thermo>::operator+=
(
    const FTMTransport<Thermo>& st
)
{

    Thermo::operator+=(st);

        // for diffusivity

        Ymd_             = st.Ymd_; 
        Xmd_             = st.Xmd_; 
        wMk_             = st.wMk_;

}


template<class Thermo>
inline void Foam::FTMTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::FTMTransport<Thermo> Foam::operator+
(
    const FTMTransport<Thermo>& st1,
    const FTMTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < small)
    {
        return FTMTransport<Thermo>
        (
            t,
            st1.coeffs_,

            // for diffusivity
            st1.Ymd_,
            st1.Xmd_,
            st1.wMk_

        );
    }
    else
    {
        //scalar Y1 = st1.Y()/t.Y();
        //scalar Y2 = st2.Y()/t.Y();

        return FTMTransport<Thermo>
        (
            t,
            st1.coeffs_,

            // for diffusivity
            st1.Ymd_,
            st1.Xmd_,
            st1.wMk_
        );
    }
}


template<class Thermo>
inline Foam::FTMTransport<Thermo> Foam::operator*
(
    const scalar s,
    const FTMTransport<Thermo>& st
)
{
    return FTMTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.coeffs_,

        // for diffusivity
        st.Ymd_,
        st.Xmd_,
        st.wMk_
    );
}


// ************************************************************************* //
