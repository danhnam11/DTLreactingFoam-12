/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Application
    chemkinToFoam

    developed by Danh Nam Nguyen and Jae Hun Lee,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description

    Utility to generate input file for using DTLreactingFoam with DTM. 
    Converts CHEMKINIII thermodynamics, transport and reaction data 
    files into OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "chemkinReader.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

Foam::string compressSpaces(const Foam::string& str)
{
    Foam::string result;
    bool lastWasSpace = false;

    forAll(str, i)
    {
        if (isspace(str[i]))
        {
            if (!lastWasSpace)
            {
                result += ' ';
                lastWasSpace = true;
            }
        }
        else
        {
            result += str[i];
            lastWasSpace = false;
        }
    }

    size_t first = result.find_first_not_of(' ');
    size_t last  = result.find_last_not_of(' ');
    return (first == Foam::string::npos) ? "" : result.substr(first, last - first + 1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(10);

    argList::validArgs.append("CHEMKIN file");
    argList::validArgs.append("CHEMKIN thermodynamics file");
    //argList::validArgs.append("dummy transport file");
    argList::validArgs.append("CHEMKIN transport file");
    argList::validArgs.append("dummy transport file");

    argList::validArgs.append("OpenFOAM chemistry file");
    argList::validArgs.append("OpenFOAM thermodynamics file");

    argList::addBoolOption
    (
        "newFormat",
        "read Chemkin thermo file in new format"
    );

    argList args(argc, argv);

    bool newFormat = args.optionFound("newFormat");

    speciesTable species;
    //chemkinReader cr(species, args[1], args[3], args[2], newFormat);
    chemkinReader cr(args[1], args[4], args[2], newFormat);
    const HashPtrTable<chemkinReader::thermoPhysics>& speciesThermo =
        cr.speciesThermo();

    dictionary thermoDict;
    thermoDict.add("species", cr.species());

// Tran.dat part
// //////////////////////////////////////////////////////////////////////

    HashTable<List<scalar>> transportMap;

    //IFstream tranStream(args[4]);
    IFstream tranStream(args[3]);
    if (!tranStream.good())
    {
        //FatalErrorInFunction << "Cannot open transport data file: " << args[4] << exit(FatalError);
        FatalErrorInFunction << "Cannot open transport data file: " << args[3] << exit(FatalError);
    }

    std::string stdLine;
    Foam::string line;
    label lineNo = 0;

    while (std::getline(tranStream().stdStream(), stdLine))
    {
        ++lineNo;

        line = stdLine;
        line.replaceAll("\t", " ");
        line.replaceAll("\r", "");

        label excl = line.find('!');
        if (excl != -1)
        {
            line = line.substr(0, excl);
        }

        if (line.empty() || line[0] == '!' || line.find_first_not_of(' ') == string::npos)
            continue;

        line = compressSpaces(line);

        std::istringstream iss(line.c_str());
        std::string species;

        label linearity;
        scalar epsilonOverKb, sigma, dipoleMoment, alpha, Zrot;

        if (!(iss >> species >> linearity >> epsilonOverKb >> sigma >> dipoleMoment >> alpha >> Zrot))
        {
            FatalErrorInFunction
                << "Failed to parse line " << lineNo << ": " << line << nl
                << "Expected format: species linearity epsilonOverKb sigma dipoleMoment alpha Zrot" << nl
                << exit(FatalError);
        }

        word speciesName(species);

        transportMap.insert(speciesName, List<scalar>{scalar(linearity), epsilonOverKb, sigma, dipoleMoment, alpha, Zrot});
    }


    // Add the species thermo formatted entries
    {
        OStringStream os;
        speciesThermo.write(os);
        dictionary speciesThermoDict(IStringStream(os.str())());

        forAll(cr.species(), i)
        {
            word speciesI = cr.species()[i];
            dictionary& speciesDict = speciesThermoDict.subDict(speciesI);

            if (!transportMap.found(speciesI))
            {
                WarningInFunction << "No transport data found for species: " << speciesI << endl;
                continue;
            }

            const List<scalar>& transportData = transportMap[speciesI];
            if (transportData.size() != 6)
            {
                FatalErrorInFunction << "Unexpected transport data size for " << speciesI << endl;
            }

            dictionary* transportDictPtr = nullptr;

            if (speciesDict.found("transport"))
            {
                transportDictPtr = const_cast<dictionary*>(&speciesDict.subDict("transport"));
            }
            else
            {
                speciesDict.add("transport", dictionary());
                transportDictPtr = &speciesDict.subDict("transport");
            }

            dictionary& transportDict = *transportDictPtr;

            transportDict.add("linearity",     transportData[0]);
            transportDict.add("epsilonOverKb", transportData[1]);
            transportDict.add("sigma",    transportData[2]);
            transportDict.add("dipoleMoment",   transportData[3]);
            transportDict.add("alpha",  transportData[4]);
            transportDict.add("Zrot", transportData[5]);

            Info << "Updated transport for " << speciesI << ": " << transportDict << endl;

        }

        thermoDict.merge(speciesThermoDict);
    }



    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    // Add elements
    forAllConstIter
    (
        HashPtrTable<chemkinReader::thermoPhysics>,
        speciesThermo,
        iter
    )
    {
        const word specieName(iter.key());

        dictionary elementsDict("elements");
        forAll(cr.specieComposition()[specieName], ei)
        {
            elementsDict.add
            (
                cr.specieComposition()[specieName][ei].name(),
                cr.specieComposition()[specieName][ei].nAtoms()
            );
        }

        thermoDict.subDict(specieName).add("elements", elementsDict);
    }

    thermoDict.write(OFstream(args[6])(), false);
    cr.reactions().write(OFstream(args[5])());

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
