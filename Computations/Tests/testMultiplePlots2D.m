clear all

% no HI, not even wall
DDFT = Test_DDFT_WallHI_Full('none');

output(1).optsPhys = DDFT.optsPhys;
output(1).optsNum  = DDFT.optsNum;
output(1).data     = DDFT.dynamicsResult;

output(1).optsPlot.lineColourDDFT = {'r'};

% just wall (i.e. position-dependent diffusion coefficient)
DDFT = Test_DDFT_WallHI_Full('noHI');

output(2).optsPhys = DDFT.optsPhys;
output(2).optsNum  = DDFT.optsNum;
output(2).data     = DDFT.dynamicsResult;

output(2).optsPlot.lineColourDDFT = {'b'};

% Full wall HI
DDFT = Test_DDFT_WallHI_Full('FullWallHI');

output(3).optsPhys = DDFT.optsPhys;
output(3).optsNum  = DDFT.optsNum;
output(3).data     = DDFT.dynamicsResult;

output(3).optsPlot.lineColourDDFT = {'g'};

% 
% % Oseen + diffusion coefficient
% DDFT = Test_DDFT_WallHI_Full('Oseen');
% 
% output(4).optsPhys = DDFT.optsPhys;
% output(4).optsNum  = DDFT.optsNum;
% output(4).data     = DDFT.dynamicsResult;
% 
% output(4).optsPlot.lineColourDDFT = {'m'};

PlotMultipleDDFTs2D(output);