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
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar Eta = 0;
    forAll(Xmd_, k)
    {
        scalar Etak = etak(wMk_[k], epsilonOverKbMk_[k], sigmaMk_[k], miuiMk_[k], p, T);        
        scalar sumPhi = 0;

        forAll(Xmd_, j)
        {
            scalar WkOj = wMk_[k]/wMk_[j];
            scalar Etaj = etak(wMk_[j], epsilonOverKbMk_[j], sigmaMk_[j], miuiMk_[j], p, T);
            scalar X = 1 + sqrt(Etak/Etaj) * pow(wMk_[j]/wMk_[k], 0.25);

            sumPhi = sumPhi + Xmd_[j]*(pow((1 + WkOj), -0.5) * pow(X, 2.0)/sqrt(8.0));
        }

        //Eta = Eta + Xmd_[k]*Etak/(sumPhi+1e-24);
        Eta = Eta + Xmd_[k]*Etak/(sumPhi);

    }
    if (Eta < 0) {Eta =  1e-16;}  
    //return Eta;  //in cgs unist
    return (0.1*Eta); // in si unit
}


template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    scalar sum1 = 0; 
    scalar sum2 = 0;
    forAll(Xmd_, k)
    {
        scalar Lamdak = 
           lamdak
           (
               wMk_[k],
               linearityMk_[k],
               epsilonOverKbMk_[k],
               sigmaMk_[k],
               miuiMk_[k],
               ZrotMk_[k],
               CpCoeffTableMk_[k],
               p,
               T
           );

           sum1 = sum1 + Xmd_[k]*Lamdak;
           sum2 = sum2 + Xmd_[k]/Lamdak;
    }

    scalar Lamda = 0.5*(sum1 + 1/(sum2 + 1e-16));
    return Lamda; //in si unit
}


// mass Diffusivity 
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::Dimix
(
    label speciei,
    const scalar p,
    const scalar T
) const
{

         scalar sum = 0.0; 
         scalar Dimix = 0.0; 

         forAll(Xmd_, j) 
         {
             if (speciei != j)
             {
             sum = sum + Xmd_[j]/
                   ( 
                       (0.00266*sqrt(pow(T,3)/(2*Mij_[j][speciei])))/
                       ((p/1.01325e5)*pow(sigmaij_[j][speciei],2)*
                        OMEGA(1, T/epsilonijOverKb_[j][speciei], deltaij_[j][speciei]) )                         
                   );
             }
         }  
         if (sum == 0.0){Dimix = 0;} 
         else {Dimix = (1-Ymd_[speciei])/(sum);}

     return Dimix/1e4;
} 


// thermal diffusivity
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::DimixT
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
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::Dij
(
    label speciei,
    label speciej,
    const scalar p,
    const scalar T
) const
{
         scalar Dij;
         if (speciei != speciej)
         {
         Dij =
             (0.00266*sqrt(pow(T,3)/(2*Mij_[speciej][speciei])))/
             ((p/1.01325e5)*pow(sigmaij_[speciej][speciei],2)*
             OMEGA(1, T/epsilonijOverKb_[speciej][speciei], deltaij_[speciej][speciei]) );
        }
        else
        {
            Dij = 1.0;
        }

     return Dij/1e4; // cgs to SI unit
}

// update coefficients for mixture calculation
template<class Thermo>
inline void Foam::preprocessingFTMTransport<Thermo>::updateTRANS
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
)
{
    Ymd_             = Ymd;
    Xmd_             = Xmd;
    epsilonijOverKb_ = epsilonijOverKb;
    deltaij_         = deltaij;
    Mij_             = Mij;
    sigmaij_         = sigmaij;

    linearityMk_     = linearityMk;
    epsilonOverKbMk_ = epsilonOverKbMk;
    sigmaMk_         = sigmaMk;
    miuiMk_          = miuiMk;
    polarMk_         = polarMk;
    ZrotMk_          = ZrotMk;
    wMk_             = wMk;
    CpCoeffTableMk_  = CpCoeffTableMk;
}


// functions to access coefficients
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::linearity() const
{   
    return linearity_;
}

template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::epsilonOverKb() const
{
    return epsilonOverKb_;
}

template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::sigma() const
{
    return sigma_;
}

template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::miui() const
{
    return miui_;
}

template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::polar() const
{
    return polar_;
}

template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::Zrot() const
{   
    return Zrot_;
}

template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::Kb() const
{
    return Kb_;
}


template<class Thermo>
inline Foam::List<Foam::scalar>
Foam::preprocessingFTMTransport<Thermo>::Ymd() const
{
    return Ymd_;
}

template<class Thermo>
inline Foam::List<Foam::scalar>
Foam::preprocessingFTMTransport<Thermo>::Xmd() const
{
    return Xmd_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::preprocessingFTMTransport<Thermo>::epsilonijOverKb() const
{
    return epsilonijOverKb_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::preprocessingFTMTransport<Thermo>::deltaij() const
{
    return deltaij_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::preprocessingFTMTransport<Thermo>::Mij() const
{
    return Mij_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::preprocessingFTMTransport<Thermo>::sigmaij() const
{
    return sigmaij_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::preprocessingFTMTransport<Thermo>::muCoeffsMk() const
{
    return muCoeffsMk_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::scalar>>
Foam::preprocessingFTMTransport<Thermo>::kappaCoeffsMk() const
{
    return kappaCoeffsMk_;
}

template<class Thermo>
inline Foam::List<Foam::List<Foam::List<Foam::scalar>>>
Foam::preprocessingFTMTransport<Thermo>::DijCoeffsMk() const
{
    return DijCoeffsMk_;
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Quadratic interpolation
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::qInterpolate
(
    scalar t,
    List<scalar> X,
    List<scalar> Y
) const
{
    scalar value = 0.0;
    scalar VAL1=0.0, VAL2=0.0, FAC1=0.0, FAC2=0.0;

    VAL1 = Y[0] + (t-X[0])*(Y[1]-Y[0]) / (X[1]-X[0]);
    VAL2 = Y[1] + (t-X[1])*(Y[2]-Y[1]) / (X[2]-X[1]);
    FAC1 = 0.5*(t-X[0]) / (X[1]-X[0]);
    FAC2 = 0.5*(X[2]-t) / (X[2]-X[1]);

    if ( t >= X[1])
    {
        value = (VAL1*FAC2+VAL2)/(1.0+FAC2);
    }
    else
    {
        value = (VAL2*FAC1+VAL1)/(1.0+FAC1);
    }
    return value;
}


//- Return collision integral from tablulated data
//- (Table IV, V Monchick and Mason, Journal of Chem. Phys. 35: 1676 (1961))
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::OMEGA
(
    int N,
    scalar TR,
    scalar DR
) const
{
     //table based on Stockmayer potentials 
     //given in Table IV, V Monchick and Mason, J. Chem. Phys. 35:1676(1961)  
    const scalar Tstar[37] = 
    { 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0, 
      1.2,  1.4,  1.6,  1.8,  2.0,  2.5,  3.0,  3.5,  4.0,  5.0,
      6.0,  7.0,  8.0,  9.0,  10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 
      25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0};  

    const scalar delta[8] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5};

    const scalar tableIV[37][8] = {    
       {4.008 , 4.002 , 4.655 , 5.52  , 6.454 , 8.214 , 9.824 , 11.31},
       {3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160},
       {2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483},
       {2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539},
       {2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918},
       {1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747},
       {1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137},
       {1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872},
       {1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657},
       {1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780},

       {1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199},
       {1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99 },
       {1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827},
       {1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7  },
       {1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592},
       {1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397},
       {0.95  , 0.9538, 0.9656, 0.9852, 1.012 , 1.082 , 1.168 , 1.265},
       {0.9131, 0.9162, 0.9256, 0.9413, 0.9626, 1.019 , 1.09  , 1.17 },
       {0.8845, 0.8871, 0.8948, 0.9076, 0.9252, 0.972 , 1.03  , 1.098},
       {0.8428, 0.8446, 0.850 , 0.859 , 0.8716, 0.9053, 0.9483, 0.9984},

       {0.813 , 0.8142, 0.8183, 0.825 , 0.8344, 0.8598, 0.8927, 0.9316},
       {0.7898, 0.791 , 0.794 , 0.7993, 0.8066, 0.8265, 0.8526, 0.8836},
       {0.7711, 0.772 , 0.7745, 0.7788, 0.7846, 0.8007, 0.822 , 0.8474},
       {0.7555, 0.7562, 0.7584, 0.7619, 0.7667, 0.78  , 0.7976, 0.8189},
       {0.7422, 0.743 , 0.7446, 0.7475, 0.7515, 0.7627, 0.7776, 0.796 },
       {0.72022,0.7206, 0.722 , 0.7241, 0.7271, 0.7354, 0.7464, 0.76  },
       {0.7025, 0.703 , 0.704 , 0.7055, 0.7078, 0.7142, 0.7228, 0.7334},
       {0.68776, 0.688, 0.6888, 0.6901, 0.6919, 0.697 , 0.704 , 0.7125},
       {0.6751, 0.6753, 0.676 , 0.677 , 0.6785, 0.6827, 0.6884, 0.6955},
       {0.664 , 0.6642, 0.6648, 0.6657, 0.6669, 0.6704, 0.6752, 0.681 },

       {0.6414, 0.6415, 0.6418, 0.6425, 0.6433, 0.6457, 0.649 , 0.653 },
       {0.6235, 0.6236, 0.6239, 0.6243, 0.6249, 0.6267, 0.629 , 0.632 },
       {0.60882,0.6089, 0.6091, 0.6094, 0.61  , 0.6112, 0.613 , 0.6154},
       {0.5964, 0.5964, 0.5966, 0.597 , 0.5972, 0.5983, 0.600 , 0.6017},
       {0.5763, 0.5763, 0.5764, 0.5766, 0.5768, 0.5775, 0.5785, 0.58  },
       {0.5415, 0.5415, 0.5416, 0.5416, 0.5418, 0.542 , 0.5424, 0.543 },
       {0.518 , 0.518 , 0.5182, 0.5184, 0.5184, 0.5185, 0.5186, 0.5187} };
     //end of Table IV
     
    const scalar tableV[37][8] = {
        {4.1   , 4.266 , 4.833 , 5.742 , 6.729 , 8.624 ,10.34  ,11.890 },
        {3.263 , 3.305 , 3.516 , 3.914 , 4.433 , 5.57  , 6.637 , 7.618 },
        {2.84  , 2.836 , 2.936 , 3.168 , 3.511 , 4.329 , 5.126 , 5.874 },
        {2.531 , 2.522 , 2.586 , 2.749 , 3.004 , 3.64  , 4.282 , 4.895 },
        {2.284 , 2.277 , 2.329 , 2.46  , 2.665 , 3.187 , 3.727 , 4.249 },
        {2.084 , 2.081 , 2.13  , 2.243 , 2.417 , 2.862 , 3.329 , 3.786 },
        {1.922 , 1.924 , 1.97  , 2.072 , 2.225 , 2.641 , 3.028 , 3.435 },
        {1.7902, 1.795 , 1.84  , 1.934 , 2.07  , 2.417 , 2.788 , 3.156 },
        {1.682 , 1.689 , 1.733 , 1.82  , 1.944 , 2.258 , 2.596 , 2.933 },
        {1.593 , 1.60  , 1.644 , 1.725 , 1.84  , 2.124 , 2.435 , 2.746 },
        
        {1.455 , 1.465 , 1.504 , 1.574 , 1.67  , 1.913 , 2.181 , 2.45  },
        {1.355 , 1.365 , 1.4   , 1.461 , 1.544 , 1.754 , 1.989 , 2.228 },
        {1.28  , 1.289 , 1.321 , 1.374 , 1.447 , 1.63  , 1.838 , 2.053 },
        {1.222 , 1.231 , 1.26  , 1.306 , 1.37  , 1.532 , 1.718 , 1.912 },
        {1.176 , 1.184 , 1.209 , 1.25  , 1.307 , 1.45  , 1.618 , 1.795 },
        {1.0933, 1.1   , 1.119 , 1.15  , 1.193 , 1.304 , 1.435 , 1.578 },
        {1.039 , 1.044 , 1.06  , 1.083 , 1.117 , 1.204 , 1.31  , 1.428 },
        {0.9996, 1.004 , 1.016 , 1.035 , 1.062 , 1.133 , 1.22  , 1.32  },
        {0.9699, 0.9732, 0.983 , 0.9991, 1.021 , 1.08  , 1.153 , 1.236 },
        {0.9268, 0.9291, 0.936 , 0.9473, 0.9628, 1.005 , 1.058 , 1.12  },
             
        {0.8962, 0.8979, 0.903 , 0.9114, 0.923 , 0.9545, 0.9955, 1.044 },
        {0.8727, 0.8741, 0.878 , 0.8845, 0.8935, 0.918 , 0.9505, 0.9893},
        {0.8538, 0.8549, 0.858 , 0.8632, 0.8703, 0.890 , 0.9164, 0.9482},
        {0.8379, 0.8388, 0.8414, 0.8456, 0.8515, 0.868 , 0.8895, 0.916 },
        {0.8243, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.89  },
        {0.8018, 0.8024, 0.8039, 0.8065, 0.810 , 0.820 , 0.8337, 0.8504},
        {0.7836, 0.784 , 0.7852, 0.7872, 0.7899, 0.7976, 0.808 , 0.8212},
        {0.7683, 0.7687, 0.7696, 0.771 , 0.7733, 0.7794, 0.788 , 0.7983},
        {0.7552, 0.7554, 0.7562, 0.7575, 0.7592, 0.764 , 0.771 , 0.7797},
        {0.7436, 0.7438, 0.7445, 0.7455, 0.747 , 0.7512, 0.757 , 0.7642},
        
        {0.71982,0.72  , 0.7204, 0.7211, 0.7221, 0.725 , 0.7289, 0.7339},
        {0.701 , 0.7011, 0.7014, 0.702 , 0.7026, 0.7047, 0.7076, 0.7112},
        {0.68545,0.6855, 0.686 , 0.686 , 0.6867, 0.6883, 0.6905, 0.693 },
        {0.6723, 0.6724, 0.6726, 0.673 , 0.6733, 0.6745, 0.676 , 0.6784},
        {0.651 , 0.651 , 0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546},
        {0.614 , 0.614 , 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147},
        {0.5887, 0.5889, 0.5894, 0.59  , 0.5903, 0.5901, 0.5895, 0.5885} };
     //end of Table V

    scalar omega = 0.0;
    if 
    ( 
        (DR < -1.0E-5) || (DR > 2.5) || (TR < 0.09) || (TR > 500) || 
        ( (DR > 1.0E-5) && (TR > 75) )
    )
    {     
        Info << " <<<<<<<<<<< ERROR: Collision integral is undefined <<<<<<<<<<" << endl;
        Info << " DR = " << DR << endl; 
        Info << " TR = " << TR << endl; 

    } 
    else
    {
        if (TR > 75)
        {
            if (N == 1)      //for OMEGA11
            {
                omega = 0.623 - 0.136E-2*TR + 0.346E-5*TR*TR -0.343E-8*TR*TR*TR;
            }
            else if (N == 2) //for OMEGA22
            {
                omega = 0.703 - 0.146E-2*TR + 0.357E-5*TR*TR -0.343E-8*TR*TR*TR; 
            }
        }
        else
        {
            int rowIndex = 0, columnIndex = 0;
            List<scalar> threePointsX(3); 
            List<scalar> threePointsY(3); 
            List<scalar> interpolatedRowData(3); 

            // Step 3: find the row index
            if (TR <= 0.2){rowIndex = 1;}
            else
            {
                for (int i =1; i < 37; i++)
                {
                    if ((TR > Tstar[i-1]) && (TR <= Tstar[i]))
                    { rowIndex = i;}
                }  
            }
            
            // Step 4
            if ( DR >= 1.0E-5)
            {
                // find the column index
                if (DR <= 0.25) {columnIndex = 1;}  
                else 
                { 
                    for (int j = 1; j < 7; j++)
                    {
                        if ((delta[j-1] < DR) && (DR <= delta[j]))
                        {columnIndex = j;} 
                    }
                }

                // do the interpolation
                // step 4.2.1. interpolate to build data for correct column of DR
                forAll(interpolatedRowData, i) 
                {
                    forAll(threePointsX, j)
                    {
                        threePointsX[j] = delta[columnIndex-1+j];
                        if (N == 1)       //for OMEGA11
                        { threePointsY[j] = tableIV[rowIndex-1+i][columnIndex-1+j];}
                        else if (N == 2)  //for OMEGA22
                        { threePointsY[j] = tableV[rowIndex-1+i][columnIndex-1+j];}                    
                    }
                    interpolatedRowData[i] = qInterpolate(DR, threePointsX, threePointsY);
                }
                // step 4.2.2. interpolate to find the data for correct row of TR
                forAll(threePointsX, i)
                {
                    threePointsX[i] = Tstar[rowIndex-1+i];
                }

                omega = qInterpolate(TR, threePointsX, interpolatedRowData);
            }
            else 
            {
                forAll(threePointsX, i)
                {
                    threePointsX[i] = Tstar[rowIndex-1+i];
                    if (N == 1)      // for OMEGA11
                    { interpolatedRowData[i] = tableIV[rowIndex-1+i][0]; }
                    else if (N == 2) // for OMEGA22
                    { interpolatedRowData[i] = tableV[rowIndex-1+i][0]; }
                }
                omega = qInterpolate(TR, threePointsX, interpolatedRowData);
            } 
        }
    }
    return omega;
}


// function return Zrot(T)
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::ZrotT
(
    scalar Zrot0,
    scalar epsilonOverKb,
    scalar T
) const
{
    // this is pi^(3/2)/2 
    const scalar c1 = 2.7841639984158539;
    // this is pi^2/4 + 2
    const scalar c2 = 4.4674011002723397;
    // this is p^(3/2)
    const scalar c3 = 5.5683279968317078;
  
    scalar D  = epsilonOverKb/T;
    scalar D0 = epsilonOverKb/298;
    scalar F = (1.0 + c1*sqrt(D0) + c2*D0 + c3*pow(D0, 1.5))/
               (1.0 + c1*sqrt(D) + c2*D + c3*pow(D, 1.5));

    return (Zrot0*F);
}


//- Temporarily calculate the dynamic viscosity of pure species
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::etak
(
    const scalar Wk,
    const scalar epsilonOverKbk,
    const scalar sigmak,
    const scalar miuik,
    const scalar p,
    const scalar T
) const
{
    scalar TR = T/epsilonOverKbk;
    scalar FDTCGS = 1.0E-18, FATCM = 1.0E+08;
    scalar Fconvert = pow(FDTCGS, 2.0)*pow(FATCM, 3.0);
    scalar DR = 0.5*Fconvert*pow(miuik, 2.0) /(1e7*Kb_* epsilonOverKbk * pow(sigmak, 3.0));
    scalar etak = 0;
    etak  = 2.6693e-5 * sqrt(Wk) * sqrt (T)/ (pow(sigmak, 2.0)*OMEGA(2, TR, DR) );
  
    return etak; // in cgs unist
}


//- Temporarily calculate the thermal conductivity [W/mK] of pure species
template<class Thermo>
inline Foam::scalar Foam::preprocessingFTMTransport<Thermo>::lamdak
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
) const
{
    scalar Pi = 3.14159265359;
    scalar Ru = 8.314510e7; //erg/(mole-K)
    scalar Rjh = 8.314510e3; //erg/(mole-K)

    scalar TR = T/epsilonOverKbk;
    scalar FDTCGS = 1.0E-18, FATCM = 1.0E+08;
    scalar Fconvert = pow(FDTCGS, 2.0)*pow(FATCM, 3.0);
    scalar DR = 0.5*Fconvert*pow(miuik, 2.0) /(1e7*Kb_* epsilonOverKbk * pow(sigmak, 3.0));  
    scalar etak = 0;
    etak  = 2.6693e-5 * sqrt(Wk) * sqrt (T)/ (pow(sigmak, 2.0)*OMEGA(2, TR, DR) );
    scalar Dkk = 0;
    Dkk = 2.6280e-3 * sqrt(pow(T, 3.0)) /
          ( p/1.01325e5 * sqrt(Wk) * pow(sigmak, 2.0) * OMEGA(1, TR, DR) );

    //scalar rho = p*Wk/(Ru*T);
    scalar PRUT = 10*p/(Ru*T);  //p in g/(cm-s^2). 1Pa = 10 g/cm-s^2
    //scalar rhoDEta = rho*Dkk/etak;
    scalar rhoDEta = PRUT*Wk*Dkk/etak;

    scalar A = 2.5 - rhoDEta;

    // calculate Cv of k-th species
    scalar Tcommon = this->Tcommon();
    List<scalar> a(7);
    if (T < Tcommon)
    {
        a = CpCoeffTablek[0];
    }
    else
    {
        a = CpCoeffTablek[1];
    }

    // Convert universal gas constant using molecular weight of species -- added by JH     
    scalar CpMCvjh = Rjh / Wk;

    // original CpMCv(p,T) function uses mixture molecular weight to retrieve R, the universal gas constant
    //scalar Cvk = ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]) - this->CpMCv(p,T);

    // Use New Cp - Cv value with species molecular weight
    scalar Cvk = ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]) - CpMCvjh;

    Cvk = 1e-3*Cvk*Wk;  //convert to mole basis, J/mol-K
    Cvk = 1e7*Cvk;      //convert to erg/(mol-K)

    //    
    scalar B = 0; 
    scalar CrOvCt = 0; // This is (Cv,rot)/(Cv,trans) ratio
    if (linearityk == 2)
    {
        B = ZrotT(Zrotk, epsilonOverKbk, T) + 2.0*(5.0/2.0 + rhoDEta)/Pi;
        CrOvCt = 1.0;
    }
    else
    {
        B = ZrotT(Zrotk, epsilonOverKbk, T) + 2.0*(5.0/3.0 + rhoDEta)/Pi;
        CrOvCt = 2.0/3.0;
    }
    
    scalar ftran = 2.5*(1.0 - 2.0*CrOvCt*A/(B*Pi));
    
    scalar frot = rhoDEta*( 1.0 + 2.0*A/(B*Pi) );
    scalar fvib = rhoDEta;
    scalar factor = 0;

    if (linearityk == 0)       // single Atom
    {
        factor = (5.0/2.0) * 1.5 * Ru;   
    }
    else if (linearityk == 1)  // linear molecule
    {
        factor = (ftran * 1.5 + frot)*Ru + fvib*(Cvk - 2.5*Ru);
    }
    else if (linearityk == 2) // non-linear molecule
    {
        factor = (ftran + frot)*1.5*Ru + fvib*(Cvk - 3.0*Ru);
    }

    //scalar lamdak = etak*factor/Wk;
    scalar lamdak = 1e-5*etak*factor/Wk; // convert to W/m-K

    return lamdak;
}


// ************************************************************************* //
