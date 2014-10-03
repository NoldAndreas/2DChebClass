clear all

% inertia
%DDFT = Test_DDFT_InertiaInfSpace_MF_osc(true);
DDFT = Test_DDFT_InertiaInfDisc_MF_osc(true);

output(1).optsPhys = DDFT.optsPhys;
output(1).optsNum  = DDFT.optsNum;
output(1).data     = DDFT.dynamicsResult;

output(1).optsPlot.lineColourDDFT = {'r'};

% diffusion
%DDFT = Test_DDFT_InertiaInfSpace_MF_osc(false);
DDFT = Test_DDFT_InertiaInfDisc_MF_osc(false);

output(2).optsPhys = DDFT.optsPhys;
output(2).optsNum  = DDFT.optsNum;
output(2).data     = DDFT.dynamicsResult;

output(2).optsPlot.lineColourDDFT = {'b'};

PlotMultipleDDFTs2D(output);