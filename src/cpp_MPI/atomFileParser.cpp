#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "atomFileParser.h"

using namespace std;

mapCollisionData collisionFileParser(string filePath, string collisionFileName, double &time, mapParticleIdToType& mapPartIdToType)
{
    //Reading Dump files
    mapCollisionData mapData;
    ifstream collisionFile;
    collisionFile.open((filePath + collisionFileName).c_str(), ifstream::in);
    if (!collisionFile.is_open())
    {
        std::cout << "Unable to open " << collisionFileName << " file" << endl;
        return mapData;
    }
    string line;
    string tmpStr;
    stringstream lineData;

    //get time value
    getline(collisionFile, line); //Read first line as ITEM: TIMESTEP
    lineData = move(stringstream(line));
    lineData >> tmpStr;
    lineData >> tmpStr;

    getline(collisionFile, line); //Read time value
    lineData = move(stringstream(line));
    lineData >> time;

    while (tmpStr.compare("ATOMS") && !collisionFile.eof())
    {
        getline(collisionFile, line);
        lineData = move(stringstream(line));
        lineData >> tmpStr;
        lineData >> tmpStr;
    }

    if (tmpStr.compare("ATOMS") || collisionFile.eof())
    {
        cout << collisionFileName << " doesn't contain ITEM:ATOMS row" << endl;
        return mapData;
    }

    while (tmpStr.compare("fz"))
        lineData >> tmpStr;

    if (tmpStr.compare("fz")) // || atomFile.eof())
    {
        cout << collisionFileName << " doesn't contain all required columns" << endl;
        return mapData;
    }

    int c_ccCount = 0;
    while (tmpStr.compare("f_fppacc"))
    {
        lineData >> tmpStr;
        c_ccCount++;
    }
    c_ccCount--;

    if (tmpStr.compare("f_fppacc")) // || atomFile.eof())
    {
        cout << collisionFileName << " doesn't contain all required columns" << endl;
        return mapData;
    }

    while (getline(collisionFile, line))
    {
        lineData = move(stringstream(line));

        //read particle id
        long int particleId = 0;
        lineData >> particleId;         
        
        //read particle type
        int particleType = 0;
        lineData >> particleType; 

        //auto mapEntry = make_pair<particleId, particleType>;
        //mapPartIdToType.insert(mapEntry);
        mapPartIdToType[particleId] = particleType;

        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore x, y & z value;
        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore ix, iy & iz value;

        collisionData cData;
        lineData >> cData.velocity[0] >> cData.velocity[1] >> cData.velocity[2]; //read vx, vy & vz value;
        lineData >> tmpStr >> tmpStr >> tmpStr;//read and ignore fx, fy & fz value;
                                          
        //read collision data
        cData.c_ccVec.resize(c_ccCount);
        for (int i = 0; i < c_ccCount; i++)
            lineData >> cData.c_ccVec[i];

        lineData >> cData.f_fpacc; //read f_fpacc value

        double radius = 0.0;
        lineData >> radius;

        auto mapIt = mapData.find(particleType);
        if (mapIt == mapData.end())
        {
            vector<collisionData> tmpVec;
            tmpVec.push_back(cData);
            auto tupleEntry = make_tuple(radius * 2, tmpVec);
            auto mapEntry = make_pair(particleType, tupleEntry);
            mapData.insert(mapEntry);
        }
        else
            get<1>(mapIt->second).push_back(cData);
    }

    //    cout << collisionFile << endl;
    //    cout << mapData.size() << endl;
    //    for(auto it = mapData.begin(); it != mapData.end(); it++)
    //        cout << (it->second).size() << endl;
    //    cout << endl;

    collisionFile.close();

    return mapData;
}

mapImpactData impactFileParser(string filePath, string impactFileName, mapParticleIdToType mapPartIdToType)
{
    //Read Collision file
    mapImpactData mapData;
    
    ifstream impactFile;
    impactFile.open((filePath + impactFileName).c_str(), ifstream::in);
    if (!impactFile.is_open())
    {
        std::cout << "Unable to open " << impactFileName << " file" << endl;
        return mapData;
    }

    string line;
    string tmpStr;
    stringstream lineData;

    getline(impactFile, line); //Read first line as ITEM: TIMESTEP
    lineData = move(stringstream(line));
    lineData >> tmpStr; //Read ITEM:
    lineData >> tmpStr; //Read TIMESTEP

    while (tmpStr.compare("ENTRIES") && !impactFile.eof())
    {
        getline(impactFile, line);
        lineData = move(stringstream(line));
        lineData >> tmpStr;
        lineData >> tmpStr;
    }

    if (tmpStr.compare("ENTRIES") || impactFile.eof())
    {
        cout << impactFileName << " doesn't contain ITEM:ENTRIES row" << endl;
        return mapData;
    }

    //int c_pwcCount = 0;
    while (tmpStr.compare("c_pwc[14]"))
    {
        lineData >> tmpStr;
        //c_pwcCount++;
    }
    //c_pwcCount--;

    if (tmpStr.compare("c_pwc[14]")) // || atomFile.eof())
    {
        cout << impactFileName << " doesn't contain all required columns" << endl;
        return mapData;
    }

    while (getline(impactFile, line))
    {
        lineData = move(stringstream(line));
        impactData iData;
        //Read velocity data: Column 1 to 6
        for(int i = 0; i < 6; i++)
        lineData >> iData.velocity[i];
        
        //Read impact type: column 7
        lineData >> iData.impactType;

        //Read & ignore mesh triangle id: column 8
        lineData >> tmpStr;

        //read particle id: column 9
        lineData >> iData.particleId; 
        
        //get particle type from id
        int particleType = 0;
        bool flag = getParticleTypeFromId(mapPartIdToType, iData.particleId, particleType);
        
        if (flag == false)
            {
                cout << "Unknown Particle Id is found in impact file; missing in collision file " << endl;
                continue;
            }
        
        //Read & ignore force value: column 10 to 12
        lineData >> tmpStr >> tmpStr >> tmpStr;

        //Read contact area: column 13
        lineData >> iData.contactArea;

        //Read overlap area: column 14
        lineData >> iData.overlapArea;

        auto mapIt = mapData.find(particleType);
        if (mapIt == mapData.end())
        {
            vector<impactData> tmpVec;
            tmpVec.push_back(iData);
            auto mapEntry = make_pair(particleType, tmpVec);
            mapData.insert(mapEntry);
        }
        else
            mapIt->second.push_back(iData);
    }

    impactFile.close();
    return mapData;
}

bool getParticleTypeFromId (mapParticleIdToType mapPartIdToType, long int particleId, int& particleType)
{
  bool retVal = false;
  auto mapIt = mapPartIdToType.find(particleId);
  if (mapIt != mapPartIdToType.end())
  {
    particleType = mapIt->second;
    retVal = true; 
  }
  return retVal;
}
