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

typedef std::tuple<double, std::vector<collisionData>> tupleDiameterAndCollisionData;
//for collision files map< particle type, tuple<vector of each row, diameter of particle>>
typedef std::map<int, tupleDiameterAndCollisionData> mapCollisionData;
//for impact files
//first value for impact with wall
//second value for impact with impeller
typedef std::pair<int, int> pairImpactData;

mapCollisionData collisionFileParser(std::string filePath, std::string collisionFileName, double &time);

pairImpactData impactFileParser(std::string filePath, std::string impactFileName);
#endif // ATOMFILEDATA_H
