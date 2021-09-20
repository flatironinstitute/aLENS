/**
 * @file ProteinConfig.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-01-04
 *
 * @copyright Copyright (c) 2019
 *
 */

#ifndef PROTEINCONFIG_HPP_
#define PROTEINCONFIG_HPP_

#include "ProteinType.hpp"

#include "KMC/lookup_table.hpp"

#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief read the proteinConfig.yaml file to read in different types of
 * proteins
 *
 * The lifetime of ProteinConfig should be the entire lifetime of the simulation
 * because types[i].LUTablePtr holds the LUT table for each protein forever
 *
 */
class ProteinConfig {
  public:
    double KBT;
    std::vector<ProteinType> types; ///< settings for different types

    std::vector<int> freeNumber; ///< free number for each type

    std::vector<std::vector<double>> fixedLocations;
    ///< fixed location for each type

    std::unordered_map<int, int> tagLookUp; ///< hash table for map and iType

    /**
     * @brief Construct a new ProteinConfig object
     * 
     * @param proteinConfigFile 
     */
    explicit ProteinConfig(std::string proteinConfigFile);

    /**
     * @brief Destroy the ProteinConfig object
     * 
     */
    ~ProteinConfig();

    /**
     * @brief display settings
     *
     */
    void echo() const;
};

#endif
