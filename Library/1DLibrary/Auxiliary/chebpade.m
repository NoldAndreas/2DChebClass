% CHEBPADE computes the poles z of an [m,n] Chebyshev-Pade approximation of a
% function with Chebyshev coefficients a. Note that a must be a column vector
% with the coefficients in the order a0, a1, a2, ...
function z = chebpade(m,n,a)
    N = length(a)-1; % number of coefficients less one
    E = chebeval([zeros(1,n);eye(N+n,n)]);
    F = repmat(chebeval([a;zeros(n,1)]),1,n);
    D = -chebpoly(E.*F);
    c = [1;D(m+2:m+n+1,:)\a(m+2:m+n+1)]; % denominator coefficients
    if (n > 1) % denominator is at least quadratic
        C = diag(0.5*ones(n-1,1),-1); % companion matrix
        C = C + diag([1;0.5*ones(n-2,1)],1);
        C(n,1:n) = -0.5*c(1:n)'/c(n+1);
        C(n,n-1) = C(n,n-1) + 0.5;
        z = eig(C); % poles
    else % denominator is linear
        z = -c(1)/c(2); % poles
    end
end