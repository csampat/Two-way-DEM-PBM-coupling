#ifndef ATOMFILEPARSER_H
#define ATOMFILEPARSER_H

#include <vector>
#include <array>
#include <map>
#include <tuple>
//#include <pair>
#include <string>

typedef struct
{
    std::array<double, 3> velocity;
    std::vector<int> c_ccVec;
    double f_fpacc;
} collisionData;

typedef struct
{
    std::array<double, 6> velocity;
    long int particleId;
    int impactType;
    double contactArea;
    double overlapArea;
} impactData;

typedef std::tuple<double, std::vector<collisionData>> tupleDiameterAndCollisionData;
//for collision files map< particle type, tuple<vector of each row, diameter of particle>>
typedef std::map<int, tupleDiameterAndCollisionData> mapCollisionData;
// for mapping of particle id with its type
typedef std::map<long int, int> mapParticleIdToType;
//for impact files
//key is particle type
//value vector of impact data
typedef std::map<int, std::vector<impactData>> mapImpactData;

bool getParticleTypeFromId (mapParticleIdToType mapPartIdToType, long int particleId, int& particleType);

mapCollisionData collisionFileParser(std::string filePath, std::string collisionFileName, double &time, mapParticleIdToType& mapPartIdToType);

mapImpactData impactFileParser(std::string filePath, std::string impactFileName, mapParticleIdToType mapPartIdToType);
#endif // ATOMFILEDATA_H
