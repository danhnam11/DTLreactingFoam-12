Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
fluidMulticomponentThermo& thermo = pThermo();
thermo.validate(args.executable(), "h");
thermo.Y();
PtrList<volScalarField>& Y = thermo.Y();

using ThermoType = preprocessingFTMTransport
        <
            species::thermo<janafThermo<perfectGas<specie> >, sensibleEnthalpy>
        >;

using MixtureType = DTMMulticomponentMixture<ThermoType>;

IOdictionary thermoDict
(
    IOobject
    (
        "physicalProperties",
        "constant",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

MixtureType mixture(thermoDict);

volScalarFieldListSlicer Yslicer(mesh, thermo.Y());
auto composition0 = Yslicer.slice(0);

ThermoType mixture_ = mixture.thermoMixture(composition0);

Foam::wordList speciesList = mixture.specieNames();
int nSpecies = mixture.specieNames().size();
const PtrList<ThermoType>& specieThermos = mixture.specieThermos();
//

Info<< "speciesList = " << speciesList << endl;
Info<< "speciesThermos = " << specieThermos << endl;
Info<< "mu of N2 = " << specieThermos[4].mu(10325, 300) << endl;

#include "createFTM.H"
