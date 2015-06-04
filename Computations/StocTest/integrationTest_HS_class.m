function integrationTest_HS_class

R = 1;

params.R = 1;

params.sigmaHS = 0.5;

HSgeom.N = [5;5]; HSgeom.L1 = 2; HSgeom.L2 = 2;
HSgeom.y2wall = 0; HSgeom.N2bound = 10; HSgeom.h = 1;
HSgeom.L2_AD = 1; HSgeom.alpha_deg = 90; HSgeom.R = R;

HS = HalfSpace_FMT(HSgeom);

Ageom.N = [10;10]; Ageom.RMin = R; Ageom.L = 2;
Ageom.shape = 'InfAnnulus';

%M_int = HS.ComputeAreaIntegrationFiniteSupport(Ageom,@RP12_2D,params,true);

M_int = HS.ComputeAreaIntegrationFiniteSupport(Ageom,@RP12_2D_noConv,params);

area = InfAnnulus(Ageom);

weights = {'HIweight_11','HIweight_12','HIweight_21','HIweight_22'};

f12      = str2func('RP12_2D');
params.HIfn = f12;
HITemp12 = HS.ComputeConvolutionFiniteSupport(area,weights,HS.Pts,params);

HITemp12(:,:,2);

M_int(~isfinite(M_int)) = 0;
M_int(:,:,1,1);

HITemp12(:,:,5) - M_int(:,:,2,2);

pause

rho = exp(-(HS.Pts.y1_kv.^2 + HS.Pts.y2_kv.^2));

val11 = M_int(:,:,1,1)*rho;
val12 = M_int(:,:,1,2)*rho;
val21 = M_int(:,:,2,1)*rho;
val22 = M_int(:,:,2,2)*rho;

[val11 val12 val21 val22 HS.Pts.y1_kv HS.Pts.y2_kv]

end