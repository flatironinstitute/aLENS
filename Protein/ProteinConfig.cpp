#include "ProteinConfig.hpp"

#include "Util/Logger.hpp"
#include "Util/YamlHelper.hpp"

ProteinConfig::ProteinConfig(std::string proteinConfigFile) {
    YAML::Node node = YAML::LoadFile(proteinConfigFile);

    readConfig(node, VARNAME(KBT), KBT, "");
    readConfig(node, VARNAME(skipBindKinetics), skipBindKinetics, "", true);

    const int nTypes = node["proteins"].size();
    types.resize(nTypes);
    freeNumber.resize(nTypes);
    fixedLocations.resize(nTypes);

    for (int i = 0; i < nTypes; i++) {
        const auto &p = node["proteins"][i];
        auto &type = types[i];

        // read type from Yaml
        type.readFromYaml(p);

        // set initial free and fixed proteins
        readConfig(p, "freeNumber", freeNumber[i], "free proteins in SimBox");
        readConfig(p, "fixedLocationPerMT", fixedLocations[i],
                   "fixed locations of proteins");

        tagLookUp[type.tag] = i;
        if (type.fixedEnd0) {
            // end0 is fixed, set probability of changing binding status to zero
            type.ko_s[0] = 0;
            type.ko_d[0] = 0;
        }
    }

    // sanity check
    for (int iType = 0; iType < nTypes; iType++) {
        if (types[iType].fixedEnd0 && freeNumber[iType] != 0) {
            spdlog::warn("Protein Type {} setting error", iType);
            spdlog::warn("Proteins with fixedEnd0 must set freeNumber = 0");
            spdlog::warn("Setting freeNumber = 0 for this type");
            freeNumber[iType] = 0;
        }
    }
}

ProteinConfig::~ProteinConfig() {
    freeNumber.clear();
    fixedLocations.clear();
    for (auto &t : types) {
        if (t.LUTablePtr)
            delete t.LUTablePtr;
    }
    types.clear();
}

void ProteinConfig::echo() const {
    printf("KBT = %g\n", KBT);
    printf("skipBindKinetics = %d\n", skipBindKinetics);
    const int nTypes = types.size();
    printf("%d Types of proteins\n", nTypes);

    for (int i = 0; i < nTypes; i++) {
        types[i].echo();
        printf("freeNumber: %d\n", freeNumber[i]);
        printf("fixedLocationPerMT: ");
        for (auto &l : fixedLocations[i]) {
            printf("%g,", l);
        }
        printf("\n");
    }
}

bool ProteinConfig::checkSkipBindKinetics() const {
    if (skipBindKinetics)
        return true;

    for (const auto &protT : types) {
        for (size_t i = 0; i < 2; i++) {
            if (protT.ko_s[i] * protT.Ka[i] * protT.eps > 0. ||
                protT.ko_s[i] > 0. ||
                protT.ko_d[i] * protT.Ke[i] * protT.eps > 0. ||
                protT.ko_d[i] > 0.)
                return false;
        }
    }

    return true;
}