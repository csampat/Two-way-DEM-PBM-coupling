#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "atomFileParser.h"

using namespace std;

mapCollisionData collisionFileParser(string filePath, string collisionFileName, double &time)
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
        //cout << collisionFileName << " doesn't contain require info" << endl;
        return mapData;
    }

    while (tmpStr.compare("fz"))
        lineData >> tmpStr;

    if (tmpStr.compare("fz")) // || atomFile.eof())
    {
        //cout << collisionFileName << " doesn't contain require info" << endl;
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
        //cout << fileName << " doesn't contain require info" << endl;
        return mapData;
    }

    while (getline(collisionFile, line))
    {
        lineData = move(stringstream(line));
        lineData >> tmpStr; //read and ignore particle id
        lineData >> tmpStr; //read particle type
        int particleType = stoi(tmpStr);
        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore x, y & z value;
        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore ix, iy & iz value;
        collisionData cData;
        lineData >> cData.velocity[0] >> cData.velocity[1] >> cData.velocity[2]; //read vx, vy & vz value;
        lineData >> tmpStr >> tmpStr >> tmpStr;                                  //read and ignore fx, fy & fz value;
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

pairImpactData impactFileParser(string filePath, string impactFileName)
{
    //Read Collision file
    pairImpactData pairData;
    pairData.first = 0;
    pairData.second = 0;
    ifstream impactFile;
    impactFile.open((filePath + impactFileName).c_str(), ifstream::in);
    if (!impactFile.is_open())
    {
        std::cout << "Unable to open " << impactFileName << " file" << endl;
        return pairData;
    }

    string line;
    string tmpStr;
    stringstream lineData;

    getline(impactFile, line); //Read first line as ITEM: TIMESTEP
    lineData = move(stringstream(line));
    lineData >> tmpStr;
    lineData >> tmpStr;

    while (tmpStr.compare("ENTRIES") && !impactFile.eof())
    {
        getline(impactFile, line);
        lineData = move(stringstream(line));
        lineData >> tmpStr;
        lineData >> tmpStr;
    }

    if (tmpStr.compare("ENTRIES") || impactFile.eof())
    {
        //cout << impactFileName << " doesn't contain require info" << endl;
        return pairData;
    }

    int impactType = 0;

    while (getline(impactFile, line))
    {
        lineData = move(stringstream(line));
        lineData >> impactType; //read impact type
        if (impactType == 0)
            pairData.first += 1;
        if (impactType == 1)
            pairData.second += 1;
    }
    impactFile.close();
    return pairData;
}
