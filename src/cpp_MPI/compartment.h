#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <vector>
#include "utility.h"
#include "parameters.h"

typedef struct
{
    arrayOfDouble2D fAll;
    arrayOfDouble2D fLiquid;
    arrayOfDouble2D fGas;
    double liquidAdditionRate;

    std::vector<double> vs;
    std::vector<double> vss;

    arrayOfDouble2D sMeshXY;
    arrayOfDouble2D ssMeshXY;

    arrayOfInt2D sAggregationCheck;
    arrayOfInt2D ssAggregationCheck;

    arrayOfInt2D sInd;
    arrayOfInt2D ssInd;

    arrayOfInt2D sIndB;
    arrayOfInt2D ssIndB;

    arrayOfDouble2D sLow;
    arrayOfDouble2D sHigh;

    arrayOfDouble2D ssLow;
    arrayOfDouble2D ssHigh;

    arrayOfInt2D sCheckB;
    arrayOfInt2D ssCheckB;

    arrayOfDouble2D diameter;

} CompartmentIn;

typedef struct
{
    arrayOfDouble2D dfAlldt;
    arrayOfDouble2D dfLiquiddt;
    arrayOfDouble2D dfGasdt;
    arrayOfDouble2D liquidBins;
    arrayOfDouble2D gasBins;
    arrayOfDouble2D internalVolumeBins;
    arrayOfDouble2D externalVolumeBins;
    arrayOfDouble4D aggregationKernel;
    arrayOfDouble4D breakageKernel;
} CompartmentOut;

typedef struct
{
    std::vector<double> DEMDiameter;
    arrayOfDouble2D DEMCollisionData;
    std::vector<double> DEMImpactData;
} CompartmentDEMIn;

typedef struct
{
    arrayOfDouble2D fAllPreviousCompartment;
    arrayOfDouble2D flPreviousCompartment;
    arrayOfDouble2D fgPreviousCompartment;
    arrayOfDouble2D fAllComingIn;
    arrayOfDouble2D fgComingIn;
} PreviousCompartmentIn;

CompartmentOut performCompartmentCalculations(PreviousCompartmentIn prevCompIn, CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double time, double timeStep);

#endif // COMPARTMENT_H
