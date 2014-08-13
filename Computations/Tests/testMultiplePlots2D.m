clear all

%DDFT = Test_DDFT_DiffusionHalfSpace_FMT_Triangle(false);

DDFT = Test_DDFT_DiffusionInfSpace_FMT(false);

output(1).optsPhys = DDFT.optsPhys;
output(1).optsNum  = DDFT.optsNum;
output(1).data     = DDFT.dynamicsResult;

output(1).optsPlot.lineColourDDFT = {'r'};

%DDFT = Test_DDFT_DiffusionHalfSpace_FMT_Triangle(true);

DDFT = Test_DDFT_DiffusionInfSpace_FMT(true);

output(2).optsPhys = DDFT.optsPhys;
output(2).optsNum  = DDFT.optsNum;
output(2).data     = DDFT.dynamicsResult;

output(2).optsPlot.lineColourDDFT = {'b'};

PlotMultipleDDFTs2D(output);