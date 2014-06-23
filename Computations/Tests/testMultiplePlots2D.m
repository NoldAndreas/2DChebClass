clear all

[data,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionInfSpace_FMT(false);

output(1) = v2struct(optsPhys,optsNum,data);

[data,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionInfSpace_FMT(true);

output(2) = v2struct(optsPhys,optsNum,data);

PlotMultipleDDFTs2D(output);