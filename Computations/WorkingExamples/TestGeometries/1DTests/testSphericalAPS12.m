clear all

[data,optsPhys,optsNum,optsPlot] = DDFT_InertiaHISphericalInfInterval();

struct(1) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = DDFT_InertiaSphericalInfInterval();

struct(2) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = DDFT_DiffusionHISphericalInfInterval();

struct(3) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = DDFT_DiffusionSphericalInfInterval();

struct(4) = v2struct(data,optsPhys,optsNum,optsPlot);


PlotMultipleDDFTs(struct);