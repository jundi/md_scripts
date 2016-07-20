/*
 * This file is modified version of a file distributed as part of the
 * GROMACS molecular simulation package.
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_DISTANCE_H
#define GMX_TRAJECTORYANALYSIS_MODULES_DISTANCE_H

#include "gromacs/trajectoryanalysis/analysismodule.h"

namespace gmx
{

namespace analysismodules
{

class DistanceInfo
{
    public:
        static const char name[];
        static const char shortDescription[];
        static TrajectoryAnalysisModulePointer create();
};

} // namespace analysismodules

} // namespace gmx

#endif
