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
#include "DP1VLU.H"

void DP1VLU
(   int L,
    int NDER, 
    double X, 
    double& YFIT, 
    List<double>& YP, 
    List<double>& A
)
{
    if (L < 0) return;
    int NDO = std::min(NDER, L);
    double VAL;
    double MAXORD = A[0] + 0.5;
    
    int K1 = static_cast<int>(MAXORD) + 1;
    int K2 = K1 + static_cast<int>(MAXORD);
    int K3 = K2 + static_cast<int>(MAXORD) + 2;
    int K4 = K3 + L + 1; 

    List<double> C(L + 1, 0.0); 

    // L IS 0
    if (L == 0) 
    {
        YFIT = A[1];
        return;
    }
    // L is 1
    if (L == 1) 
    {
        C[1] = A[2];
        YFIT = A[1] + (X - A[2]) * C[1];
        if (NDER >= 1) 
        {
            YP[1] = C[1];
        }
        return;
    }

    // L IS GREATER THAN 1
    int NDP1 = NDO + 1;
    int K3P1 = K3 + 1;
    int K4P1 = K4 + 1;
    int LP1 = L + 1;
    int LM1 = L - 1;
    int ILO = K3 + 3;
    int IUP = K4 + NDP1;

    for (int i = ILO-1; i <= IUP-1; ++i)
    {
        A[i] = 0.0;
    }

    double DIF = X - A[LP1-1];
    int KC = K2 + LP1;
    A[K4P1-1] = A[KC-1];
    A[K3P1-1] = A[KC - 1-1] + DIF * A[K4P1-1];
    A[K3 + 2-1] = A[K4P1-1];

    // EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
    for (int I = 1; I <= LM1; ++I)
    {
        int IN = L - I-1;
        int INP1 = IN + 1;
        int K1I = K1 + INP1;
        int IC = K2 + IN;

        double DIF = X - A[INP1];
        VAL = A[IC] + DIF * A[K3P1-1] - A[K1I] * A[K4P1-1];

        if (NDO <= 0) break;

        for (int N = 1; N <= NDO; ++N) 
        {
            int K3PN = K3P1 + N-1;
            int K4PN = K4P1 + N-1;
            YP[N-1] = DIF * A[K3PN] + N * A[K3PN - 1] - A[K1I] * A[K4PN];
        }

        // SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
	for (int N = 1; N <= NDO; ++N)
	{
            int K3PN = K3P1 + N-1;
            int K4PN = K4P1 + N-1;
	    A[K4PN] = A[K3PN];
	    A[K3PN] = YP[N-1];
	}

	A[K4P1-1] = A[K3P1-1];
	A[K3P1-1] = VAL;

    }
    YFIT = VAL;

    Info << "DP1VLU function call was successful." << endl;
}

