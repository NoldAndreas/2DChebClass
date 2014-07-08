clear all

geom.N = [5;8]; geom.L1 = 5; geom.L2 = 10;

aBox = Box(geom);

aBox.ComputeIndices;

M = aBox.M;
Ind = aBox.Ind;

Ind.normalFinite;

nnz(Ind.finite1)
size(Ind.normalFinite1)
nnz(Ind.finite2)
size(Ind.normalFinite2)