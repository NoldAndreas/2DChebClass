function A = barychebevalMatrix(x,xx)

% BARYCHEBEVAL computes the values ff taken by a rational interpolant with barycentric
% weights w, interpolation points x, and function values f at the points xx using
% formula (3.1). Note that x must be arranged in increasing order.

% DNS: Assume (3.1) is a reference to Andreas' own writeup.
% DNS: x are the usual Chebyshev points
% DNS: xx are the 

% This seems to follow Barycentric interpolation. Tref has a paper:
% Barycentric Lagrange Interpolation, Jean-Paul Berrut Lloyd N. Trefethen
% SIAM REVIEW Vol. 46, No. 3, pp. 501–517
 
n = length(x);      % DNS: set n as the number of discretisation points
nn = length(xx);    % DNS: set nn as the number of new points interpolated onto
        
%Insert from Nikos  
w = (-1).^(0:(n-1))'.*[0.5;ones(n-2,1);0.5];
% DNS: Vector is 1/2 -1 1 ... -1 1 -1/2 (this is 1./c in Tref book)
% DNS: NOTE: These are NOT ClenCurt weights. This is the same w as in my
% Cheb Diff Matrices function
       
[sorted,indices] = sort([x(2:n);xx(1:nn)]);
% DNS: This sorts the vectors x and xx into ascending order, the indices
% are kept track of in indices, so that sorted = [x(2:n);xx(1:nn)](indices)
        
%indices<=n-1 labels entries belonging to x 
%(i.e. having an index <=n-1) with 1, otherwise 0

%cumsum(indices<=n-1)+1
%adds elements up           % DNS: cumsum keeps vector structure, but
                            % changes the elements to be the cumulative sum
                        
sorted(indices) = cumsum(indices<=n-1)+1;
 
%to better understand        
%[sorted,[x(2:n);xx(1:nn)]]
        
%array with the indices of the neighboring left points in x
% for each point in xx 
index = sorted(n:n+nn-1);
        
maska = (x(index)==xx); % mark points in xx = x
maskb = (maska==0);     % mark points in xx ~= x
xxx   = xx(maskb);      % store in xxx points s.t. xx ~= x
        
% numer = zeros(length(xxx),1);   %array with length of unknown points
% DNS: Was commented out by Andreas
denom = zeros(length(xxx),1);    %   "
A     = zeros(nn,n);
        
if(~isempty(xxx))
    for i=1:n    %go through all points in x
        temp = w(i)./(xxx-x(i));        % DNS: Paper p507: temp = c(j)./xdiff
        A(maskb,i) = temp;
        denom = denom + temp;
    end
end
% DNS: This loop is like in the paper mentioned above, p506-507
% c = [1/2; ones(n-1,1); 1/2].*(-1).^((0:n)'); --- which is basically the
% same as w given above, but with n=n-1 in c. Then xx is defined, then
% numer = zeros(size(xx)); --- similar
% denom = zeros(size(xx)); --- similar
% for j = 1:n+1; xdiff = xx-x(j); temp = c(j)./xdiff --- Same essentially
% numer = numer + temp*f(j); --- this line isn't there
% denom = demon + temp; end; --- this is

A(maskb,:)  =  diag(1./denom) *  A(maskb,:) ;

for i=1:length(maska)
    if(maska(i)==1)
        A(i,index(i)) = 1;
    end
end
        

% DNS: These were commented out by Andreas - bit like p510 of paper
%        ff = zeros(nn,1);
%        ff(maskb) = numer./denom;
%        ff(maska) = f(index(maska));

end
