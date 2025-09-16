forAll(Y,i)
{
    Info << " =========== print speciesName ==============  " << endl;
    Info << " ith = " << i << ";--> speciesName =  " << speciesList[i] << endl;

    // Fit Viscosity
    DPOLFT(N, TLog, muLog[i], W, MAXDEG, NDEG, EPS, R, IERR, muA);
    if (IERR == 1)
    {
        Info << "Successfully computed the polynomial for ." << endl;
    }
    else
    {
        Info << "Error occurred. IERR = " << IERR << endl;
    }
    DPCOEF(MAXDEG, C, muA, muTC[i]);

    Info << "Taylor series coefficients for mu:" << endl;

    forAll(muTC[i],index)
    {
        Info << "muTC(" << index << ") = " << muTC[i][index] << endl;
    }
    Info << endl;

    //Fit Thermal Conducvity
    DPOLFT(N, TLog, kappaLog[i], W, MAXDEG, NDEG, EPS, R, IERR, kappaA);
    if (IERR == 1)
    {
        Info << "Successfully computed the polynomial for ." << endl;
    }
    else
    {
        Info << "Error occurred. IERR = " << IERR << endl;
    }
    DPCOEF(MAXDEG, C, kappaA, kappaTC[i]);

    Info << "Taylor series coefficients for kappa:" << endl;

    forAll(kappaTC[i],index)
    {
        Info << "kappaTC(" << index << ") = " << kappaTC[i][index] << endl;
    }
    Info << endl;

    //Fit Binary Diffusion Coefficients
    forAll(speciesList, j)
    {
        DPOLFT(N, TLog, DijLog[i][j], W, MAXDEG, NDEG, EPS, R, IERR, DijA[j]);
        if (IERR == 1)
        {
            Info << "Successfully computed the polynomial for ." << endl;
        }
        else
        {
            Info << "Error occurred. IERR = " << IERR << endl;
        }
        DPCOEF(MAXDEG, C, DijA[j], DijTC[i][j]);

        Info << "Taylor series coefficients for Dij:" << endl;

        forAll(DijTC[i][j],index)
        {
            Info << "DijTC(" << index << ") = " << DijTC[i][j][index] << endl;
        }
        Info << endl;
    }

    Info << endl;

    propertyReader.printTranCoeff(runTime, speciesList[i], muTC[i], kappaTC[i], DijTC[i]);

}
