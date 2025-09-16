const std::string inputFile = "constant/thermo.DTM";
const std::string outputFile = "constant/thermo.FTM";

std::ifstream in(inputFile);
std::ofstream out(outputFile);

if (!in.is_open())
{
    std::cerr << "Failed to open " << inputFile << std::endl;
    return 1;
}

HashTable<label> speciesToIndex;
forAll(speciesList, i)
{
    speciesToIndex.insert(speciesList[i], i);
}

std::string line;
std::string trimmed;
std::string currentSpecies;
std::string pendingSpecies;
bool waitingForOpenBrace = false;
bool insideElements = false;
bool insideTransport = false;
bool writingTransport = false;

while (std::getline(in, line))
{
    trimmed = line;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

    if (!waitingForOpenBrace &&
        trimmed != "" &&
        trimmed != "{" &&
        trimmed != "}" &&
        trimmed != "elements" &&
        trimmed != "specie" &&
        trimmed != "thermodynamics" &&
        trimmed != "transport" &&
        trimmed.find(' ') == std::string::npos)
    {
        pendingSpecies = trimmed;
        waitingForOpenBrace = true;
    }

    if (waitingForOpenBrace && trimmed == "{")
    {
        currentSpecies = pendingSpecies;
        waitingForOpenBrace = false;
    }

    if (trimmed == "transport")
    {
        insideTransport = true;
        writingTransport = true;
        continue;
    }

    if (insideTransport && trimmed == "}")
    {
        insideTransport = false;
        writingTransport = false;
        continue;
    }

    if (writingTransport)
    {
        continue;
    }

    if (trimmed == "elements")
    {
        insideElements = true;
    }

    if (insideElements && trimmed == "}")
    {
        insideElements = false;
        out << line << std::endl;

        if (speciesToIndex.found(currentSpecies))
        {
            label speciesIndex = speciesToIndex[currentSpecies];
            Info << "✅ inserting transport for species: " << currentSpecies << endl;
            
	    out << "    transport" << std::endl;
            out << "    {" << std::endl;

            // muCoeffs
            std::ostringstream mu;
            mu << "        muCoeffs        (";
            for (int i = 0; i < muTC[speciesIndex].size(); ++i)
            {
                mu << muTC[speciesIndex][i];
                if (i < muTC[speciesIndex].size() - 1) mu << " ";
            }
            mu << ");";
            out << mu.str() << std::endl;

            // kappaCoeffs
            std::ostringstream kappa;
            kappa << "        kappaCoeffs     (";
            for (int i = 0; i < kappaTC[speciesIndex].size(); ++i)
            {
                kappa << kappaTC[speciesIndex][i];
                if (i < kappaTC[speciesIndex].size() - 1) kappa << " ";
            }
            kappa << ");";
            out << kappa.str() << std::endl;

            // DijCoeffs
            out << "        DijCoeffs       (" << std::endl;
            for (int j = 0; j < speciesList.size(); ++j)
            {
                std::ostringstream dij;
                dij << "                            (";
                for (int d = 0; d < DijTC[speciesIndex][j].size(); ++d)
                {
                    dij << DijTC[speciesIndex][j][d];
                    if (d < DijTC[speciesIndex][j].size() - 1) dij << " ";
                }
                dij << ")";
                out << dij.str() << std::endl;
            }
            out << "        );" << std::endl;
            out << "    }" << std::endl;
        }
        else
        {
            Info << "⚠️  Skipping unknown species: " << currentSpecies << endl;
        }
        continue;
    }

    out << line << std::endl;
}

in.close();
out.close();
