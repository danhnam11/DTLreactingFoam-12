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
#include "PropertyReader.H"

PropertyReader::PropertyReader
(
    const Foam::fvMesh& mesh,
    const Foam::wordList& speciesList
)
:
    mesh(mesh),
    speciesNames(speciesList)
{

}



void PropertyReader::generateTranCoeff
(
    const Foam::wordList& speciesList    
)
{
    std::string dataDirectory = "constant/";

    std::string fileName = dataDirectory + "TranportCoeffDict";
    std::ofstream outFile(fileName);

    if (!outFile.is_open())
    {
        FatalErrorInFunction << "Error opening output file!" << exit(FatalError);
    }

    // Write OpenFoam Header content
    outFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  6\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "\\*---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    version     2.0;\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    location    \"constant\";\n";
    outFile << "    object      TranCoeff;\n";
    outFile << "}\n";
    outFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
    outFile << "\n";
    // End of header

    outFile << "\n";
    outFile << "species\n";
    outFile << "(\n";
    forAll(speciesList, i)
    {
	outFile << "    "<<  speciesList[i] <<"\n";
    }
    outFile << ");\n";

    outFile.close();
}

void PropertyReader::printTranCoeff
(
    const Foam::Time& runTime,
    const Foam::word& species,
    List<scalar>& muTCi,
    List<scalar>& kappaTCi,
    List<List<scalar>>& DijTCi

)
{

    std::string dataDirectory = "constant/";

    std::string fileName = dataDirectory + "TranportCoeffDict";
    std::ofstream outFile(fileName, std::ios::app); // Append mode

    outFile << "\n";

    outFile << species << "\n";
    outFile << "{\n";
    outFile << "    " << "fittingTranCoeff\n";
    outFile << "    {\n";

    // Print muCoeffs in the desired format
    outFile << "        " << "muCoeffs    (";
    forAll(muTCi, degree)
    {
        outFile << muTCi[degree];
        if (degree < muTCi.size() - 1) // Add a space if it's not the last element
        {
            outFile << " ";
        }
    }

    outFile << ");\n"; // Close the coefficients with a semicolon

    // Print kappaCoeffs in the desired format
    outFile << "        " << "kappaCoeffs    (";
    forAll(kappaTCi, degree)
    {
        outFile << kappaTCi[degree];
        if (degree < kappaTCi.size() - 1) // Add a space if it's not the last element
        {
            outFile << " ";
        }
    }

    // Print DijCoeffs in the desired format
    outFile << ");\n"; // Close the coefficients with a semicolon

    outFile << "        DijCoeffs    (\n";

    forAll(DijTCi, j)
    {
        outFile << "                        (";

        forAll(DijTCi[j], degree)
        {

            outFile << DijTCi[j][degree];

            if (degree < DijTCi[j].size() - 1)
            {
                outFile << " ";
            }
        }

        outFile << ")\n";
    }

    outFile << "                     );\n";


    outFile << "    }\n";
    outFile << "};\n";

}
