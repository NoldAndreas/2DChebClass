% CHEBPOLY uses the fast Fourier transform to find M+1 Chebyshev coefficients
% A of a function, given its values U at M+1 Chebyshev points. The operation
% is applied to each column of U.
function A = chebpoly(U)
    M = size(U,1)-1;
    U(M+2:2*M,:) = U(M:-1:2,:);
    A = real(ifft(U)); A = A(1:M+1,:); A(2:M,:) = 2*A(2:M,:);
end