#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <regex>
#include <map>

main ()
{
    unsigned int nNonNeutrons(0), nMuons(1), nElectrons(0), nProtons(1), nPiPlus(0), nPiMinus(0), nPhotons(3), nKaonPlus(0), nKaonMinus(0);
    int nuNuanceCode(5010);

    //logic to deconstruct new nuance codes >= 5000 (if there are kaons we ignore this bit: kaons are not included in the interactionType enum anyway)
    if (nuNuanceCode >= 5000 && nKaonPlus == 0 && nKaonMinus == 0)
    {   
        //deconstruct custom nuance code into CC/NC and mode
        int mutableNuNuanceCode(nuNuanceCode);
        std::vector<int> baseTenDeconstruction;

        while (mutableNuNuanceCode > 0)
        {   
            baseTenDeconstruction.push_back(mutableNuNuanceCode % 10);
            mutableNuNuanceCode /= 10; 
        }   

        int CCNC(baseTenDeconstruction.at(2));
        int mode(10 * baseTenDeconstruction.at(1) + baseTenDeconstruction.at(0));
        std::string interactionTypeString("");

        //define CC/NC
        if (CCNC == 1)
            interactionTypeString.append("CC");
        else
            interactionTypeString.append("NC");

        std::cout << "mode: " << mode << std::endl;

        std::map<int, std::string> modeToInteractionTypeMap = {{0, "QEL_"}, {1, "RES_"}, {2, "DIS_"}, {3, "COH_"}, {10, "MEC_"}};
        interactionTypeString.append(modeToInteractionTypeMap.at(mode));


        std::vector<std::pair<int, std::string>> particleCountsWithNames = {{nMuons, "MU_"}, {nElectrons, "E_"}, {nProtons, "P_"}, 
                                                                            {nPiPlus, "PIPLUS_"}, {nPiMinus, "PIMINUS_"}, {nPhotons, "PHOTON_"}};

        for (const auto &pair : particleCountsWithNames)
        {   
            for (int i = 0; i < pair.first; ++i)
                interactionTypeString.append(pair.second);
        }   

        //remove trailing underscore
        interactionTypeString.erase(interactionTypeString.size() - 1); 

        //map two photons to PIZERO, as before
        //std::string doublePhotonString("PHOTON_PHOTON");
        //index = interactionTypeString.find(doublePhotonString, index);
        //interactionTypeString = interactionTypeString.replace(index, doublePhotonString.length(), "PIZERO");

        interactionTypeString = std::regex_replace(interactionTypeString, std::regex("PHOTON_PHOTON"), "PIZERO");

        std::cout << interactionTypeString << std::endl; 
    }   
}
