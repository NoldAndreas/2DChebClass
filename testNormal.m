clear all


disp('Box');
geom.N = [5;8]; geom.L1 = 5; geom.L2 = 10;

aBox = Box(geom);
aBox.ComputeIndices;

Ind = aBox.Ind;

x = aBox.Pts.y1_kv;
y = aBox.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]

disp('InfSpace');
geom.N = [5;8]; geom.L1 = 5; geom.L2 = 10;

aInfSpace = InfSpace(geom);
aInfSpace.ComputeIndices;

Ind = aInfSpace.Ind;

x = aInfSpace.Pts.y1_kv;
y = aInfSpace.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]

disp('HalfSpace');
PhysArea = struct('N',[20,20],'L1',2,'L2',2,'y2wall',0.,...
                          'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);
theta_CS     = PhysArea.alpha_deg*pi/180;
optsHS       = PhysArea;
optsHS.alpha = theta_CS;
aHalfSpace   = HalfSpace_FMT(optsHS,0.5);

aHalfSpace.ComputeIndices;

Ind = aHalfSpace.Ind;

x = aHalfSpace.Pts.y1_kv;
y = aHalfSpace.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]

disp('Disc');
geom.N = [5;8]; geom.R = 5;

aDisc = Disc(geom);
aDisc.ComputeIndices;

Ind = aDisc.Ind;

x = aDisc.Pts.y1_kv;
y = aDisc.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]

disp('InfDisc');
geom.N = [5;8]; geom.L = 2;

aInfDisc = InfDisc(geom);
aInfDisc.ComputeIndices;

Ind = aInfDisc.Ind;

x = aInfDisc.Pts.y1_kv;
y = aInfDisc.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]

disp('InfAnnulus');
geom.N = [5;8]; geom.L = 2; geom.RMin = 2;

aInfAnnulus = InfAnnulus(geom);
aInfAnnulus.ComputeIndices;

Ind = aInfAnnulus.Ind;

x = aInfAnnulus.Pts.y1_kv;
y = aInfAnnulus.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]

disp('Annulus');
geom.N = [5;8]; geom.RMin = 2; geom.RMax = 5;

aAnnulus = Annulus(geom);
aAnnulus.ComputeIndices;

Ind = aAnnulus.Ind;

x = aAnnulus.Pts.y1_kv;
y = aAnnulus.Pts.y2_kv;

[Ind.normalFinite1*x Ind.normalFinite1*y]
[Ind.normalFinite2*x Ind.normalFinite2*y]
