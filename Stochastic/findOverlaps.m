function [pairs,nij,Rij] = findOverlaps(x,sigma,dim)

%overlaps = zeros(length(x)/dim);

[Rij,nij] = getRij(x,x,dim);

% add in diagonal elements to avoid counting self-overlaps
Rij = Rij+2*diag(diag(sigma));

minSep = sigma;

sepTest = Rij - minSep;

overlaps = zeros(size(Rij));
overlaps(sepTest<0) = 1;

[row,col] = find(triu(overlaps));

pairs = [row,col];

end

