#include <cmath>

#include "kernel.h"
#include "utility.h"
#include "liggghtsData.h"
using namespace std;

#define DUMP(varName) dumpData(varName, #varName)
#define DUMP2D(varName) dump2DData(varName, #varName)
#define DUMP3D(varName) dump3DData(varName, #varName)

#define DUMPCSV(varName) dumpCSV(varName, #varName)
#define DUMP2DCSV(varName) dump2DCSV(varName, #varName)
#define DUMP4DCSV(varName) dump4DCSV(varName, #varName)

#define PI 3.14159

arrayOfDouble4D DEMDependentAggregationKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, arrayOfDouble2D externalLiquidContent, double timeStep)
{
    arrayOfDouble4D aggregationKernel = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //Initialization
    arrayOfDouble4D collisionFrequency = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
    arrayOfDouble4D collisionEfficiency = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //Collision Frequency (from 2D Number of Collisions)
   
    arrayOfDouble2D fAll = compartmentIn.fAll;

    arrayOfDouble2D DEMCollisionData = compartmentDEMIn.DEMCollisionData;

    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    collisionFrequency[s1][ss1][s2][ss2] = (DEMCollisionData[s2][ss2] * timeStep) / TIMESTEPDEM;

    // Constant Collision Efficiency
    // FROM Barrasso, Tamrakar, Ramachandran. Procedia Engineering 102 (2015) 1295�1304. (p. 1298)
    //    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
    //        for (int ss1  = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
    //            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
    //                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
    //                    collisionEfficiency[s1][ss1][s2][ss2] = COLLISIONEFFICIENCYCONSTANT;

    // FROM Sen, Barrasso, Singh, Ramachandran. Processes 2014, 2, 89-111. (p. 96)
    // double collisionEfficiencyConstant=0.01;
    double criticalExternalLiquid = 0.2;
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                {
                    bool flag1 = (fAll[s1][ss1] >= 0.0) && (fAll[s2][ss2] >= 0.0);
                    bool flag2 = (externalLiquidContent[s1][ss1] >= criticalExternalLiquid) && (externalLiquidContent[s2][ss2] >= criticalExternalLiquid);
                    if (flag1 && flag2)                        
                        collisionEfficiency[s1][ss1][s2][ss2] = COLLISIONEFFICIENCYCONSTANT;
                }
    
    //Aggregation Kernel Calculation
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    aggregationKernel[s1][ss1][s2][ss2] = AGGREGATIONKERNELCONSTANT * collisionFrequency[s1][ss1][s2][ss2] * collisionEfficiency[s1][ss1][s2][ss2];

    return aggregationKernel;
}

arrayOfDouble4D DEMDependentBreakageKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double timeStep)
{
    arrayOfDouble4D breakageKernel = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //INITIALIZATION
    //arrayOfDouble2D impactFrequency = getArrayOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    arrayOfDouble2D fAll = compartmentIn.fAll;
    vector<double> numberOfImpacts = compartmentDEMIn.DEMImpactData;
    vector<double> impactFrequency = compartmentDEMIn.DEMImpactData;   
    double Ubreak = 0.0;
    //Impact Frequency (from 1D Number of Impacts)
    for (int s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
        for (int ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            for (int i = 0; i < NUMBEROFDEMBINS - 1; i++)
            {
                if (fAll[s][ss] > 0.0)
                    //impactFrequency[s][ss] = (numberOfImpacts[i] * timeStep) / TIMESTEPDEM;
                    impactFrequency[i] = (numberOfImpacts[i] * timeStep) / TIMESTEPDEM;
            }

    // Critical velocity for breakage
    Ubreak = (2 * CRITICALSTOKESDEFNUMBER / SOLIDDENSITY) * (9 / 8) * (pow((1 - INITIALPOROSITY),2) / pow(INITIALPOROSITY,2)) * (9 / 16) * (BINDERVISCOSITY / compartmentIn.diameter[0][0]); 
    
    liggghtsData* lData = liggghtsData::getInstance();
    vector<double> velocity = lData->getFinalDEMVelocity();
    if ((velocity).size() == 0)
    {
        cout << "Velocity data is missing in LIGGGHTS impact file" << endl;
        return breakageKernel;
    }
    
    double sum = 0.0;
    int size1 = velocity.size();

    for (int i = 0; i < size1; i++)
    	sum += velocity[i];

    double averageVelocity = sum / NUMBEROFDEMBINS;
    double stdDevVelocity = 0.0;
    double varianceVelocity = 0.0;
    for (int i = 0; i < size1; ++i)
    {
        varianceVelocity += pow((velocity[i] - averageVelocity), 2) / 10;
    }
    
    stdDevVelocity = sqrt(varianceVelocity);
    double intVelocity = 0.0;

    vector<double> probablityOfVelocity;
    for (int i = 0; i < size1; i++)
    {
     	probablityOfVelocity[i] = (1 / (velocity[i] * sqrt(2 * PI) * stdDevVelocity)) * exp(-((log(velocity[i]) - averageVelocity) / (2 * varianceVelocity)));
     	if(velocity[i] > Ubreak)
     		intVelocity += probablityOfVelocity[i] * velocity[i];
    }

    //DUMP(impactFrequency);
    
    //Breakage Kernel Calculation
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    breakageKernel[s1][ss1][s2][ss2] = impactFrequency[ss1] * intVelocity * BREAKAGEKERNELCONSTANT;
                    //breakageKernel[s1][ss1][s2][ss2] = impactFrequency[ss1] * BREAKAGEPROBABILITY * BREAKAGEKERNELCONSTANT;
                    //breakageKernel[s1][ss1][s2][ss2] = impactFrequency[s1][ss1] * BREAKAGEPROBABILITY;

    return breakageKernel;
    
}
