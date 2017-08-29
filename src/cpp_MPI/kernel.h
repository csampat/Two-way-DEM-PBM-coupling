#ifndef KERNEL_H
#define KERNEL_H

#include "utility.h"
#include "parameters.h"
#include "compartment.h"

arrayOfDouble4D DEMDependentAggregationKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, arrayOfDouble2D externalLiquidContent, double timeStep);
arrayOfDouble4D DEMDependentBreakageKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double timeStep);

#endif // #define KERNEL_H
