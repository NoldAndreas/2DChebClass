clear all

[data,optsPhys,optsNum,optsPlot] = DDFT_InertiaInfInterval();

struct(1) = v2struct(data,optsPhys,optsNum,optsPlot);

[data,optsPhys,optsNum,optsPlot] = NS_InfInterval();

struct(2) = v2struct(data,optsPhys,optsNum,optsPlot);

PlotMultipleDDFTs(struct);