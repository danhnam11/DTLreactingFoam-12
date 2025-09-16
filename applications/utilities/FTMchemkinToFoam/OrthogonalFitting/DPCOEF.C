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
#include "DPCOEF.H"
#include "DP1VLU.H"

void DPCOEF
(
    int L, double C,
    List<double>& A,
    List<double>& TC
)
{
    int LL = std::abs(L);
    int LLP1 = LL + 1;

    double YFIT;
    List<double> YP(std::max(0, L), 0.0);
    DP1VLU(LL, LL, C, YFIT, YP, A);
 
    TC[0] = YFIT;

    for (int I = 1; I <= LL; ++I)
    {
        TC[I] = YP[I-1];
    }

    if (LL < 2) 
    {
        if (L >= 0) 
        {
	    return;
        } 
        else 
        {
        }
        return;
    }

    double FAC = 1.0;
    for (int I = 3; I <= LLP1; ++I)
    {
        FAC *= (I - 1);
        TC[I - 1] /= FAC;
    }
/*
    int NR = LLP1 / 2;
    for (int I = 0; I < NR; ++I)
    {
        double SAVE = TC[I];
        int NEW = LLP1 - I - 1;
        TC[I] = TC[NEW];
        TC[NEW] = SAVE;
    }
*/
}
