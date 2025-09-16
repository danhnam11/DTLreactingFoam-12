Info<< "speciesList = " << speciesList << endl;
Info<< "speciesThermos = " << specieThermos << endl;
Info<< "mu of N2 = " << specieThermos[4].mu(10325, 300) << endl;

Info<< "Reading thermophysical properties\n" << endl;


IOdictionary FittingDict
(
    IOobject
    (
        "FittingDict",
        "constant",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

int numberOfMesh = FittingDict.lookupOrDefault<int>("numberOfMesh", 101);
scalar TMin        = FittingDict.lookupOrDefault<scalar>("minT", 300.0);
scalar TMax        = FittingDict.lookupOrDefault<scalar>("maxT", 3000.0);
scalar pressure    = FittingDict.lookupOrDefault<scalar>("pressure", 101325.0);

int N = numberOfMesh;
List<scalar> W(N, 1.0);
int MAXDEG = 3;
int NDEG;
double C = 0.0;
double EPS = 0.0;
int IERR;

// Generate the variables for fitting

List<scalar> R(N, 0.0);
List<scalar> muA(3 * N + 3 * MAXDEG + 3, 0.0);
List<scalar> kappaA(3 * N + 3 * MAXDEG + 3, 0.0);
List<List<scalar>> DijA(nSpecies);
forAll(DijA, i)
{
    DijA[i] = List<scalar>(3 * N + 3 * MAXDEG + 3, 0.0);
}

List<List<scalar>> muTC(nSpecies, List<scalar>(MAXDEG + 1, 0.0));
List<List<scalar>> kappaTC(nSpecies, List<scalar>(MAXDEG + 1, 0.0));
List<List<List<scalar>>> DijTC(nSpecies);
forAll(DijTC, i)
{
    DijTC[i].resize(nSpecies);
    forAll(DijTC[i], j)
    {
        DijTC[i][j].resize(MAXDEG + 1, 0.0);
    }
}

List<List<scalar>> muList(nSpecies, List<scalar>(N, 0.0));
List<List<scalar>> kappaList(nSpecies, List<scalar>(N, 0.0));
List<List<List<scalar>>> DijList(nSpecies);
forAll(DijList, i)
{
    DijList[i].resize(nSpecies);

    forAll(DijList[i], j)
    {
        DijList[i][j].resize(N, 0.0);
    }
}

// Generate Log List
List<List<scalar>> muLog = muList;
List<List<scalar>> kappaLog = kappaList;
List<List<List<scalar>>> DijLog = DijList;

List<scalar> TList(N);
List<scalar> pList(N);

List<scalar> TLog = TList;

forAll(TList, dataI)
{
    TList[dataI] = TMin + (TMax - TMin) * scalar(dataI) / scalar(N - 1);
    pList[dataI] = pressure;
}

// Assign mu & kappa of specie i to muList & kappaList 

forAll(TList, dataI)
{
    TLog[dataI] = Foam::log(TList[dataI]);
}

forAll(specieThermos, i)
{
    forAll(muList[i], celli)
    {
        scalar muVal     = specieThermos[i].mu(pList[celli], TList[celli]);
        scalar kappaVal  = specieThermos[i].kappa(pList[celli], TList[celli]);

        muList[i][celli]     = muVal;
        kappaList[i][celli]  = kappaVal;

        muLog[i][celli]      = Foam::log(muVal);
        kappaLog[i][celli]   = Foam::log(kappaVal);

    }
}

// Assign Dij of specie i & j to DijList
forAll(TList, celli)
{

    forAll(DijList, i)
    {
        forAll(DijList[i], j)
        {
            scalar DijVal = mixture_.Dij(i, j, pList[celli], TList[celli]);
            DijList[i][j][celli] = DijVal;
            DijLog[i][j][celli]  = Foam::log(DijVal);

        }
    }

}

PropertyReader propertyReader(mesh, speciesList);
propertyReader.generateTranCoeff(speciesList);

