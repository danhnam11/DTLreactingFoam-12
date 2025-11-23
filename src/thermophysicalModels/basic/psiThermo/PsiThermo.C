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

#include "PsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// calculate thermo variables using original transport models in OF
template<class BaseThermo>
void Foam::PsiThermo<BaseThermo>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);
/*
    Info<< "celli: " << celli << " [THERMO mixture]" << nl
        << "  muCoeffsMk: "    << thermoMixture.muCoeffsMk()    << nl
        << "  kappaCoeffsMk: " << thermoMixture.kappaCoeffsMk() << nl
        << "  DijCoeffsMk: "   << thermoMixture.DijCoeffsMk()   << nl
        << endl;
*/
        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture =
            this->transportMixture(composition, thermoMixture);

        TCells[celli] = thermoMixture.The
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                phe[facei] = thermoMixture.he(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                pT[facei] = thermoMixture.The(phe[facei], pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
    }
}


// calculate thermo variables using DTM
// only use for constructor since it is relatively expensive
template<class BaseThermo>
void Foam::PsiThermo<BaseThermo>::initialize()
{
    Info << "[Expensive!] using Detailed Transport Model WITHOUT coTHERM " << endl;        
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = this->Wmix_.primitiveFieldRef();

    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture = this->transportMixture(composition, thermoMixture);
/*
    if (celli % 100 == 0)   // 100개 셀마다 출력
    {
        Info<< "celli: " << celli << nl
            << "  p:     " << pCells[celli] << nl
            << "  T:     " << TCells[celli] << nl
            << "  mu:    " << muCells[celli] << nl
            << "  kappa: " << kappaCells[celli] << nl
            << "  muCoeffsMk: "    << transportMixture.muCoeffsMk()    << nl
            << "  kappaCoeffsMk: " << transportMixture.kappaCoeffsMk() << nl
            << "  DijCoeffsMk: "   << transportMixture.DijCoeffsMk()   << nl
            << endl;
    }
*/
        TCells[celli] = thermoMixture.The
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        WmixCells[celli]  = thermoMixture.W();
        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);

        //Info << "JH TEST 0: before mu & kappa call " << endl;

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);

        //Info << "JH TEST 0-1: after mu call " << endl;
        kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);

        //Info << "JH TEST 1: passed mu & kappa call " << endl;

        forAll(this->Dimix_, i)
        {

            //Info << "JH TEST 1-1: before Dimix call " << endl;

            this->Dimix_[i].primitiveFieldRef()[celli] 
          = transportMixture.Dimix(i, pCells[celli], TCells[celli]);

            //Info << "JH TEST 1-2: after Dimix call " << endl;

            this->DimixT_[i].primitiveFieldRef()[celli]
          = transportMixture.DimixT(i, pCells[celli], TCells[celli]);

            //Info << "JH TEST 1-3: after DimixT call " << endl;

            this->heList_[i].primitiveFieldRef()[celli]
          = specieThermos_[i].he(pCells[celli], TCells[celli]); 

            //Info << "JH TEST 1-4: after heList call " << endl;
        }
        //Info << "JH TEST 2: passed mu & kappa call " << endl;

    }

    volScalarField::Boundary& pBf =this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = this->Wmix_.boundaryFieldRef();

    //Info << "JH TEST 3: passed celli loop" << endl;

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];  

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                phe[facei] = thermoMixture.he(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);

                forAll(this->Dimix_, i)
                {
                    this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    this->heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].he(pp[facei], pT[facei]);
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                pT[facei] = thermoMixture.The(phe[facei], pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();   
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
           
                forAll(this->Dimix_, i)
                {
                    this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    this->heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].he(pp[facei], pT[facei]);
                }
            }
        }
    }
}



// calculate thermo variables using DTM with coTHERM algorithm 
// This can significantly reduce computational time
template<class BaseThermo>
void Foam::PsiThermo<BaseThermo>::calculateUsingCoTHERM()
{
    Info << "[Note!] using Detailed Transport Models + CoTHERM" << endl;
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = this->Wmix_.primitiveFieldRef();

    // For CoTHERMCount
    scalarField& CountCells = this->coTHERMStepCount_.primitiveFieldRef();

    // old fields
    const scalarField& TCellsOld = this->T_.oldTime();
    const scalarField& PCellsOld = this->p_.oldTime();

    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture = this->transportMixture(composition, thermoMixture);

        TCells[celli] = thermoMixture.The
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        this->residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli];
        scalar coDeltaT = mag(TCells[celli] - TCellsOld[celli]);
        scalar coDeltaP = mag(pCells[celli] - PCellsOld[celli]);

        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        WmixCells[celli]  = thermoMixture.W();

        const bool exceedCount = (CountCells[celli] >= this->maxCoTHERMStepCount_);
  
        if 
        ( 
            (coDeltaT <= this->epsilonT_) && 
            (this->flagSpecies_[celli] < 1.0 ) &&
            (coDeltaP <= this->epsilonP_) &&
            (!exceedCount)             
        )
        {
            // If the temperature is unchanged, don't need to calculate 
            // thermophysical properties, just copy from old time fields 
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 1.0;
            // start counter     
            CountCells[celli] += 1.0;

            CpCells[celli] = this->Cp_.oldTime()[celli];
            CvCells[celli] = this->Cv_.oldTime()[celli];
            muCells[celli] = this->mu_.oldTime()[celli];
            kappaCells[celli] = this->kappa_.oldTime()[celli];
    
            forAll(this->Dimix_, i)
            {
                this->Dimix_[i].primitiveFieldRef()[celli] 
              = this->Dimix_[i].oldTime()[celli];
    
                this->DimixT_[i].primitiveFieldRef()[celli]
              = this->DimixT_[i].oldTime()[celli];
    
                this->heList_[i].primitiveFieldRef()[celli]
              = this->heList_[i].oldTime()[celli];
            }
        }
        else if 
        (
            (coDeltaT <= this->epsilonT_) && 
            (this->flagSpecies_[celli] < 1.0 ) &&
            !(coDeltaP <= this->epsilonP_) &&
            (!exceedCount)               
        )
        {
            // only re-calculate Dimix
            // for other thermophysical properties, copy from old time
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 0.5;     
            // reset counter
            CountCells[celli] = 0.0; 

            CpCells[celli] = this->Cp_.oldTime()[celli];
            CvCells[celli] = this->Cv_.oldTime()[celli];
            muCells[celli] = this->mu_.oldTime()[celli];
            kappaCells[celli] = this->kappa_.oldTime()[celli];
    
            forAll(this->Dimix_, i)
            {
                this->Dimix_[i].primitiveFieldRef()[celli] 
              = transportMixture.Dimix(i, pCells[celli], TCells[celli]);
    
                this->DimixT_[i].primitiveFieldRef()[celli]
              = this->DimixT_[i].oldTime()[celli];
    
                this->heList_[i].primitiveFieldRef()[celli]
              = this->heList_[i].oldTime()[celli];
            }
        }
        else
        {
            // otherwise, calculate all thermophysical properties again 
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 0.0;
            // reset counter
            CountCells[celli] = 0.0; 

            CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
            CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
            muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
            kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);
      
            forAll(this->Dimix_, i)
            {
                this->Dimix_[i].primitiveFieldRef()[celli] 
              = transportMixture.Dimix(i, pCells[celli], TCells[celli]);
    
                this->DimixT_[i].primitiveFieldRef()[celli]
             = transportMixture.DimixT(i, pCells[celli], TCells[celli]);
    
                this->heList_[i].primitiveFieldRef()[celli]
              = specieThermos_[i].he(pCells[celli], TCells[celli]); 
            }
        }
    }

    volScalarField::Boundary& pBf =this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = this->Wmix_.boundaryFieldRef(); 
    volScalarField::Boundary& CountBf  = this->coTHERMStepCount_.boundaryFieldRef();

    // old fields
    const volScalarField::Boundary& TBfOld = this->T_.oldTime().boundaryField();
    const volScalarField::Boundary& pBfOld = this->p_.oldTime().boundaryField();
    const volScalarField::Boundary& CpBfOld = this->Cp_.oldTime().boundaryField();
    const volScalarField::Boundary& CvBfOld = this->Cv_.oldTime().boundaryField();
    const volScalarField::Boundary& muBfOld = this->mu_.oldTime().boundaryField();    
    const volScalarField::Boundary& kappaBfOld = this->kappa_.oldTime().boundaryField();


    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];    
        fvPatchScalarField& pCount  = CountBf[patchi];

        // old fields
        const fvPatchScalarField& pTOld = TBfOld[patchi]; 
        const fvPatchScalarField& ppOld = pBfOld[patchi];
        const fvPatchScalarField& pCpOld = CpBfOld[patchi];
        const fvPatchScalarField& pCvOld = CvBfOld[patchi];        
        const fvPatchScalarField& pmuOld = muBfOld[patchi];        
        const fvPatchScalarField& pkappaOld = kappaBfOld[patchi];
        //

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                scalar coDeltaTp = mag(pT[facei] - pTOld[facei]);
                scalar coDeltaPp = mag(pp[facei] - ppOld[facei]);                

                phe[facei] = thermoMixture.he(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();         

                const bool exceedCountP = (pCount[facei] >= this->maxCoTHERMStepCount_);

                if 
                ( 
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    (coDeltaPp <= this->epsilonP_)  &&
                    (!exceedCountP)                                      
                )
                {
                    // don't need to re-calculate thermophysical properties
                    // just copy from old time fields
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 1.0;
                    // start counter
                    pCount[facei] += 1.0;

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = this->Dimix_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = this->DimixT_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryField()[patchi][facei];
                    }
                }
                else if 
                (
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    !(coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                      
                )
                {
                    // only re-calculate Dimix
                    // for other thermophysical properties, copy from old time
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.5;
                    // reset counter
                    pCount[facei] = 0.0;

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = this->DimixT_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryField()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate all thermophysical properties again 
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.0;
                    // reset counter
                    pCount[facei] = 0.0;

                    pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                    pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                    pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                    pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.DimixT(i, pp[facei], pT[facei]);
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].he(pp[facei], pT[facei]);
                    }
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                scalar coDeltaTp = mag(pT[facei] - pTOld[facei]);
                scalar coDeltaPp = mag(pp[facei] - ppOld[facei]);                
                
                pT[facei] = thermoMixture.The(phe[facei], pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W(); 

                const bool exceedCountP = (pCount[facei] >= this->maxCoTHERMStepCount_);

                if 
                ( 
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    (coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                                         
                )
                {
                    // don't need to re-calculate thermophysical properties
                    // just copy from old time fields
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 1.0;
                    // start counter
                    pCount[facei] += 1.0;

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = this->Dimix_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = this->DimixT_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryField()[patchi][facei];
                    }
                }
                else if 
                (
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    !(coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                          
                )
                {
                    // only re-calculate Dimix
                    // for other thermophysical properties, copy from old time
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.5;
                    // reset counter
                    pCount[facei] = 0.0;

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = this->DimixT_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryField()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate all thermophysical properties again     
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.0;
                    // reset counter
                    pCount[facei] = 0.0;

                    pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                    pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                    pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                    pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                    = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                    = transportMixture.DimixT(i, pp[facei], pT[facei]);
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                    = specieThermos_[i].he(pp[facei], pT[facei]);
                    }
                }
            }
        }
    }
}



// calculate thermo variables using DTM with coTHERM algorithm, 
// but only using T criterion, only internally in our Lab
// This also can significantly reduce computational time
template<class BaseThermo>
void Foam::PsiThermo<BaseThermo>::calculateUsingCoTHERMOnlyT()
{
    Info << "[Note!] using Detailed Transport Model + CoTHERM only T " << endl;    
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = this->Wmix_.primitiveFieldRef();

    // old fields
    const scalarField& TCellsOld = this->T_.oldTime();

    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture = this->transportMixture(composition, thermoMixture);

        TCells[celli] = thermoMixture.The
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        this->residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli];

        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        WmixCells[celli]  = thermoMixture.W();

        if ( mag(TCells[celli] - TCellsOld[celli] ) < this->epsilonOnlyT_ )   
        {
            // If the temperature is unchanged, don't need to calculate 
            // thermophysical properties, just copy from old time fields 
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 1.0;      

            CpCells[celli] = this->Cp_.oldTime()[celli];
            CvCells[celli] = this->Cv_.oldTime()[celli];
            muCells[celli] = this->mu_.oldTime()[celli];
            kappaCells[celli] = this->kappa_.oldTime()[celli];
    
            forAll(this->Dimix_, i)
            {
                this->Dimix_[i].primitiveFieldRef()[celli] 
              = this->Dimix_[i].oldTime()[celli];
    
                this->DimixT_[i].primitiveFieldRef()[celli]
              = this->DimixT_[i].oldTime()[celli];
    
                this->heList_[i].primitiveFieldRef()[celli]
              = this->heList_[i].oldTime()[celli];
            }
        }
        else
        {
            // otherwise, calculate thermophysical properties again
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 0.0;    

			CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
			CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
			muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
			kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);
      
			forAll(this->Dimix_, i)
			{
				this->Dimix_[i].primitiveFieldRef()[celli] 
			  = transportMixture.Dimix(i, pCells[celli], TCells[celli]);
	
				this->DimixT_[i].primitiveFieldRef()[celli]
			  = transportMixture.DimixT(i, pCells[celli], TCells[celli]);
	
				this->heList_[i].primitiveFieldRef()[celli]
			  = specieThermos_[i].he(pCells[celli], TCells[celli]); 
			}
        }
    }

    volScalarField::Boundary& pBf =this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = this->Wmix_.boundaryFieldRef(); 

    // old fields
    const volScalarField::Boundary& TBfOld = this->T_.oldTime().boundaryField();
    const volScalarField::Boundary& CpBfOld = this->Cp_.oldTime().boundaryField();
    const volScalarField::Boundary& CvBfOld = this->Cv_.oldTime().boundaryField();
    const volScalarField::Boundary& muBfOld = this->mu_.oldTime().boundaryField();    
    const volScalarField::Boundary& kappaBfOld = this->kappa_.oldTime().boundaryField();


    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];    

        // old fields
        const fvPatchScalarField& pTOld = TBfOld[patchi]; 
        const fvPatchScalarField& pCpOld = CpBfOld[patchi];
        const fvPatchScalarField& pCvOld = CvBfOld[patchi];        
        const fvPatchScalarField& pmuOld = muBfOld[patchi];        
        const fvPatchScalarField& pkappaOld = kappaBfOld[patchi];
        //


        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];

				phe[facei] = thermoMixture.he(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();         

                if ( mag(pT[facei] - pTOld[facei] ) < this->epsilonOnlyT_ )
                {
                    // If the temperature is unchanged, don't need to calculate 
                    // thermophysical properties, just copy from old time fields 
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 1.0;      

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = this->Dimix_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = this->DimixT_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryField()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate thermophysical properties again     
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.0;      

					pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
					pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
					pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
					pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
	
					forAll(this->Dimix_, i)
					{
						this->Dimix_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.Dimix(i, pp[facei], pT[facei]);
	
						this->DimixT_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.DimixT(i, pp[facei], pT[facei]);
	
						this->heList_[i].boundaryFieldRef()[patchi][facei]
					  = specieThermos_[i].he(pp[facei], pT[facei]);
					}
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                
				pT[facei] = thermoMixture.The(phe[facei], pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W(); 

                if ( mag(pT[facei] - pTOld[facei] ) < this->epsilonOnlyT_ )
                {
                    // If the temperature is unchanged, don't need to calculate 
                    // thermophysical properties, just copy from old time fields 
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 1.0;       

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(this->Dimix_, i)
                    {
                        this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = this->Dimix_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = this->DimixT_[i].oldTime().boundaryField()[patchi][facei];
    
                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryField()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate thermophysical properties again 
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.0;       
                       
					pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
					pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
					pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
					pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
	
					forAll(this->Dimix_, i)
					{
						this->Dimix_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.Dimix(i, pp[facei], pT[facei]);
	
						this->DimixT_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.DimixT(i, pp[facei], pT[facei]);
	
						this->heList_[i].boundaryFieldRef()[patchi][facei]
					  = specieThermos_[i].he(pp[facei], pT[facei]);
					}
                }
            }
        }
    }
}



// calculate thermo variables using DTM + generate binary diffusivity coeffs
// only use this function in preprocessing steps for FTM model
template<class BaseThermo>
void Foam::PsiThermo<BaseThermo>::calculateTransportPreProcessing()
{
    Info << "[Expensive!] using Detailed Transport Model for preprocessing FTM " << endl;    
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = this->Wmix_.primitiveFieldRef();

    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture = this->transportMixture(composition, thermoMixture);

        TCells[celli] = thermoMixture.The
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        WmixCells[celli]  = thermoMixture.W();
        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);

        forAll(this->Dimix_, i)
        {
            this->Dimix_[i].primitiveFieldRef()[celli] 
          = transportMixture.Dimix(i, pCells[celli], TCells[celli]);

            this->DimixT_[i].primitiveFieldRef()[celli]
          = transportMixture.DimixT(i, pCells[celli], TCells[celli]);

            this->heList_[i].primitiveFieldRef()[celli]
          = specieThermos_[i].he(pCells[celli], TCells[celli]); 

        }

        // for binary diffusion coefficients
        forAll(this->Dij_, i)
        {
            forAll(this->Dij_[i],j)
            {
                this->Dij_[i][j].primitiveFieldRef()[celli]
              = transportMixture.Dij(i, j, pCells[celli], TCells[celli]);
            }
        }
    }

    volScalarField::Boundary& pBf =this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = this->Wmix_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];  

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                phe[facei] = thermoMixture.he(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);

                forAll(this->Dimix_, i)
                {
                    this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    this->heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].he(pp[facei], pT[facei]);
                }

                // for binary diffusion coefficients
                forAll(this->Dij_, i)
                {
                    forAll(this->Dij_[i], j)
                    {
                        this->Dij_[i][j].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dij(i, j, pp[facei], pT[facei]);
                    }
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                pT[facei] = thermoMixture.The(phe[facei], pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();   
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
           
                forAll(this->Dimix_, i)
                {
                    this->Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    this->DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    this->heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].he(pp[facei], pT[facei]);
                }

                // for binary diffusion coefficients
                forAll(this->Dij_, i)
                {
                    forAll(this->Dij_[i], j)
                    {
                        this->Dij_[i][j].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dij(i, j, pp[facei], pT[facei]);
                    }
                }

            }
        }
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::PsiThermo<BaseThermo>::PsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName),
    specieThermos_(this->specieThermos()), //
    mode_(0) //
{

    // check and setup coTHERM mode - Nam
    if 
    (
        this->usingDetailedTransportModel_ && 
        !this->usingCoTHERM_ && 
        !this->usingCoTHERMOnlyT_ &&
        !this->usingPreProcessingFTM_       
    )
    {
        //modeName_ = "DetailedModels";        
        mode_ = 1;
    }
    else if 
    (
        this->usingDetailedTransportModel_ && 
        this->usingCoTHERM_ && 
        !this->usingCoTHERMOnlyT_ &&
        !this->usingPreProcessingFTM_       
    )
    {
        //modeName_ = "coTHERM";        
        mode_ = 2;
    }
    else if
    (
        this->usingDetailedTransportModel_ && 
        !this->usingCoTHERM_ && 
        !this->usingCoTHERMOnlyT_ &&
        this->usingPreProcessingFTM_        
    )
    {
        //modeName_ = "preprocessingCoTHERM";
        mode_ = 3;
    }
    else if 
    (
        this->usingDetailedTransportModel_ && 
        !this->usingCoTHERM_ && 
        this->usingCoTHERMOnlyT_ && 
        !this->usingPreProcessingFTM_       
    )
    {
        //modeName_ = "CoTHERMonlyT";
        mode_ = 4;
    }
    else
    {
        //modeName_ = "originalOpenFOAM";
        mode_ = 0;
    }
    // Nam 

    // calculate();  // original 
    Info << "[This is in constructor] thermophysical properties initialization! "  << endl;    
    initialize();    // Nam for DTM
    Info << "[This is in constructor] end of thermophysical properties initialization! "  << endl;

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::PsiThermo<BaseThermo>::~PsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::PsiThermo<BaseThermo>::correct()
{
    if (BaseThermo::debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    // calculate(); // original 

    // Nam
    switch(mode_)
    {
        case 1 : 
            initialize();            
            break; 
            
        case 2 : 
            calculateUsingCoTHERM();
            break;      
            
        case 3 : 
            calculateTransportPreProcessing();
            break; 

        case 4 : 
            calculateUsingCoTHERMOnlyT();
            break; 
            
        default:
            calculate();
            break; 
    }

    if (BaseThermo::debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
