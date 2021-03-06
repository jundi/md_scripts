/*
 * This file is modified version of a file distributed as part of the
 * GROMACS molecular simulation package.
 */

#include "distance.h"

#include <string>

#include <gromacs/trajectoryanalysis.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

class Distance : public TrajectoryAnalysisModule
{
    public:
        Distance();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        Selection	                         sel_;
        Selection                                ref_;
        std::string                              fnXYZ_;
        std::string                              fnZ_;
        std::string                              fnAbsZ_;
        AnalysisData                             xyz_;
        AnalysisData                             z_;
        AnalysisData                             absz_;
};

Distance::Distance()
    : TrajectoryAnalysisModule(DistanceInfo::name, DistanceInfo::shortDescription)
{
    registerAnalysisDataset(&xyz_, "xyz");
    registerAnalysisDataset(&z_, "z");
    registerAnalysisDataset(&absz_, "absz");
}


void
Distance::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
	"[THISMODULE] calculates distances between a position and a set of",
	"positions as a function of time.",
    };

    options->setDescription(desc);

    options->addOption(FileNameOption("oxyz").filetype(eftPlot).outputFile()
			   .store(&fnXYZ_).defaultBasename("distxyz")
			   .description("Distance components as function of time"));
    options->addOption(FileNameOption("oz").filetype(eftPlot).outputFile()
                           .store(&fnZ_).defaultBasename("distz")
                           .description("Distance z-component as function of time"));
    options->addOption(FileNameOption("oabsz").filetype(eftPlot).outputFile()
                           .store(&fnAbsZ_).defaultBasename("absdistz")
                           .description("Absolute distance z-component as function of time"));
    options->addOption(SelectionOption("select").store(&sel_).required()
                           .description("Group to calculate distances for"));
    options->addOption(SelectionOption("ref").store(&ref_).required()
                           .description("Reference position"));
}


/*! \brief
 * Checks that selections conform to the expectations of the tool.
 */
void checkReference(const Selection &ref)
{
   if (ref.posCount() != 1)
   {
       std::string message = formatString(
                   "Selection '%s' does not evaluate into one position."
                   "(there are %d positions)",
                   ref.name(), ref.posCount());
       GMX_THROW(InconsistentInputError(message));
   }
}


void
Distance::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation         & /*top*/)
{
    const int distCount = sel_.posCount();
    xyz_.setColumnCount(0, distCount * 3);
    z_.setColumnCount(0, distCount);
    absz_.setColumnCount(0, distCount);

    if (!fnXYZ_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnXYZ_);
        plotm->setTitle("Distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        // TODO: Add legends? (there can be a massive amount of columns)
        xyz_.addModule(plotm);
    }

    if (!fnZ_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnZ_);
        plotm->setTitle("Distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        // TODO: Add legends? (there can be a massive amount of columns)
        z_.addModule(plotm);
    }

    if (!fnAbsZ_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAbsZ_);
        plotm->setTitle("Distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        // TODO: Add legends? (there can be a massive amount of columns)
        absz_.addModule(plotm);
    }
}


void
Distance::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                       TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle   xyzHandle  = pdata->dataHandle(xyz_);
    AnalysisDataHandle   zHandle  = pdata->dataHandle(z_);
    AnalysisDataHandle   abszHandle  = pdata->dataHandle(absz_);

    checkReference(ref_);

    xyzHandle.startFrame(frnr, fr.time);
    zHandle.startFrame(frnr, fr.time);
    abszHandle.startFrame(frnr, fr.time);
    xyzHandle.selectDataSet(0);
    zHandle.selectDataSet(0);
    abszHandle.selectDataSet(0);
    for (int i = 0; i < sel_.posCount(); i++)
    {
        const SelectionPosition &p1 = ref_.position(0);
        const SelectionPosition &p2 = sel_.position(i);
        rvec                     dx;
        if (pbc != NULL)
        {
            pbc_dx(pbc, p2.x(), p1.x(), dx);
        }
        else
        {
            rvec_sub(p2.x(), p1.x(), dx);
        }
        bool bPresent = p1.selected() && p2.selected();

        xyzHandle.setPoints(i*3, 3, dx);
	zHandle.setPoint(i, dx[2], bPresent);
	for(unsigned int i = 0; i < 3; i++) {
	  if(dx[i] < 0)dx[i] *= -1;
	}
        abszHandle.setPoint(i, dx[2], bPresent);
    }
    xyzHandle.finishFrame();
    zHandle.finishFrame();
    abszHandle.finishFrame();
}


void
Distance::finishAnalysis(int /*nframes*/)
{
}


void
Distance::writeOutput()
{
    printf("Number of positions in selection:  %d\n", sel_.posCount());
}

}       // namespace

const char DistanceInfo::name[]             = "distance";
const char DistanceInfo::shortDescription[] =
    "Calculate distances between pairs of positions";

TrajectoryAnalysisModulePointer DistanceInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Distance);
}

} // namespace analysismodules

} // namespace gmx

int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<gmx::analysismodules::Distance>(argc, argv);
}

