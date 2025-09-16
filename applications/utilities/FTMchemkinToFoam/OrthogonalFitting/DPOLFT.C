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
#include "DPOLFT.H"

void DPOLFT
(
    int N,
    const List<double>& X,
    const List<double>& Y,
    const List<double>& W,
    int MAXDEG,
    int& NDEG,
    double EPS,
    List<double>& R,
    int& IERR,
    List<double>& A
)
{
    int M = N;
    if (M == 0 || MAXDEG < 0 || M < (MAXDEG + 1))
    {
        IERR = 2; // Invalid parameters
        return;
    }

    A[0] = MAXDEG;
    
    int K1 = MAXDEG + 1;
    int K2 = K1 + MAXDEG;
    int K3 = K2 + MAXDEG + 2;
    int K4 = K3 + M;
    int K5 = K4 + M;

    double W11 = 0.0;
    double W1 = 0.0;

    double ETST = 0.1;

    // Initialized A array and W11
    for (int i = 0; i < M; ++i)
    {
        int K4PI = K4 + i; // for Fortran index
        A[K4PI] = 1.0;
        W11 += W[i];
    }

    // Compute the fit of degree zero
    double TEMD1 = 0.0;
    for (int i = 0; i < M; ++i)
    {
        int K4PI = K4 + i;
        TEMD1 += W[i] * Y[i] * A[K4PI];
    }
    TEMD1 /= W11;
    A[K2] = TEMD1;

    double SIGJ = 0.0;
    for (int i = 0; i < M; ++i)
    {
        int K4PI = K4 + i;
        R[i] = TEMD1 * A[K4PI];
        double diff = Y[i] - R[i];
        SIGJ += W[i] * diff * diff;
    }

    int J = 0;
    while (J <= MAXDEG)
    {
        ++J;
        int JP1 = J + 1;
        int K1PJ = K1 + J;
        int K2PJ = K2 + J;

        // Compute new B coefficient Except when J == 1
        if (J > 1)
        {
            A[K1PJ-1] = W11 / W1;
        }
        TEMD1 = 0.0;

        // Compute new A coefficient
        for (int i = 0; i < M; ++i)
        {
            int K4PI = K4 + i;
            TEMD1 += X[i] * W[i] * A[K4PI] * A[K4PI];
	}
        A[JP1-1] = TEMD1 / W11;

        // EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
        W1 = W11;
        W11 = 0.0;
        for (int i = 0; i < M; ++i)
        {
            int K3PI = K3 + i;
            int K4PI = K4 + i;
            double TEMP = A[K3PI];
            A[K3PI] = A[K4PI];
            A[K4PI] = (X[i] - A[JP1-1]) * A[K3PI] - A[K1PJ-1] * TEMP;
            W11 += W[i] * A[K4PI] * A[K4PI];
        }

        // GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE PRECISION
        TEMD1 = 0.0;
        for (int i = 0; i < M; ++i)
        {
	    //Info << "TEMD1 is:" << TEMD1 << endl;
            int K4PI = K4 + i;
            int K5PI = K5 + i;
            double TEMD2 = W[i] * ((Y[i] - R[i]) - A[K5PI]) * A[K4PI];
            TEMD1 += TEMD2;
        }
        TEMD1 /= W11;
        A[K2PJ] = TEMD1;

        // UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND ACCUMULATE SUM OF SQUARES
        SIGJ = 0.0;
        for (int i = 0; i < M; ++i)
        {
            int K4PI = K4 + i;
            int K5PI = K5 + i;
            double TEMD2 = R[i] + A[K5PI] + TEMD1 * A[K4PI];
            R[i] = TEMD2;
            A[K5PI] = TEMD2 - R[i];
            SIGJ += W[i] * std::pow(Y[i] - R[i] - A[K5PI], 2);
        }

        if (EPS == 0.0) // Handling the case when EPS is exactly 0
        {
            //Info << " Control Statements 26 has been used" << endl;
            if (MAXDEG == J)
            {
                IERR = 1; // Set error code for max degree reached
                NDEG = J; // Update the degree
                A[K3-1] = NDEG; // Update the A array with the current degree

 		Info << "DPOLFT function call was successful." << endl;
                return; // Exit the function
            }

            continue;
        }

        // CHECK IF RMS ERROR CRITERION IS SATISFIED (INPUT EPS > 0)
        else if (EPS > 0.0)
        {
            Info << " Control Statements 27 has been used" << endl;
            if (SIGJ <= ETST)
            {
                IERR = 1; // Set error code
                NDEG = J; // Update the degree
                A[K3-1] = NDEG; // Update the A array with the current degree
                return; // Exit the function
            }

            // CHECK IF MAXIMUM DEGREE HAS BEEN REACHED
            if (MAXDEG == J)
            {
                IERR = 3; // Set error code for max degree reached
                NDEG = MAXDEG; // Set the degree to max degree
                A[K3-1] = NDEG; // Update the A array with the current degree
                return; // Exit the function
            }
        }
    }
    Info << "DPOLFT function call was successful." << endl;
}

