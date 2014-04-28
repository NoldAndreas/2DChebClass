clear all

data =  TestFMT_DDFT_DiffusionInfSpace_HI_switch(0);
data_HI = TestFMT_DDFT_DiffusionInfSpace_HI_switch(1);

struct(1) = data;
struct(2) = data_HI;

PlotMultipleDDFTs2D(struct);