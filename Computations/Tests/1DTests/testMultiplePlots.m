clear all

[data,optsPhys,optsNum,optsPlot] = DDFT_InertiaInfInterval();

struct(1) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = DDFT_DiffusionInfInterval();

struct(2) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = DDFT_DiffusionHIInfInterval();

struct(3) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = DDFT_InertiaHIInfInterval();

struct(4) = v2struct(data,optsPhys,optsNum,optsPlot);

PlotMultipleDDFTs(struct);