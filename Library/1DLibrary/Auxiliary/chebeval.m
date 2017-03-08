% CHEBEVAL uses the fast Fourier transform to find the values U of a function
% at N+1 Chebyshev points, given its N+1 Chebyshev coefficients A. The
% operation is applied to each column of A.
function U = chebeval(A)
    N = size(A,1)-1;
    A(2:N,:) = A(2:N,:)/2; A(N+2:2*N,:) = A(N:-1:2,:);
    U = real(fft(A)); U = U(1:N+1,:);
end