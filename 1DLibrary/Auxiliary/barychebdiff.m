function Diff = barychebdiff(x,order)
% BARYCHEBDIFF constructs the first, second, third, and fourth order differentiation
% matrices D1, D2, D3, and D4 associated with a rational interpolant with
% barycentric weights w and interpolation points x using formula (3.28).

% DNS: This function is as in my IFM code to compute the Chebyshev
% differentiation matrices. Surprisingly, this is thus based on Nikos' code
% rather than Trefethan's (which Andreas seems to use more often). Also
% oddly, to use Trefethan's cheb.m in my IFM code, I needed to make the change
% N->N-1 in all code after defining N = length(x); (and x if appropriate).

% DNS: I don't know what formula (3.28) is in the comment above.
% (Probably just chebyshev points like my code
% x = sin(pi*(N-1:-2:1-N)'/(2*(N-1)));
% ).
% In fact, reading the files using this, x is defined in ClenCurt, and x
% there is defined as
% theta = pi*(0:N)'/N;
% x = cos(theta);
% which is actually as in Tref (cheb.m)
% DNS: Note however, that after x is defined like this, in the main code
% Andreas uses xR = flipdim(xR,1); to swap the points to increasing, i.e. [-1,1]

% DNS: Barycentric - point located at the centre of mass
% DNS: Basically the same as the chebdiff in my IFM code, but with w already specified.
% DNS: NOTE: This w is NOT the ClenCurt weights. In Andreas' main code, it
% is of the form (-1).^(0:NR)'.*[0.5;ones(NR-1,1);0.5]; whereas in my
% droplet code it is of the form w = (-1).^(0:n-1)'.*[0.5;ones(n-2,1);0.5];

% I have decided to put w back in here rather than have two seperate instances of it
% in the main code for r and theta

    Diff = struct();
    if(nargin == 1)
        order = 4;
    end    

    n = length(x);      % DNS: set n as the number of discretisation points
    
    if(n==1)
        Diff.Dx    = 0;
        Diff.DDx   = 0;
        Diff.DDDx  = 0;
        Diff.DDDDx = 0;
        return;
    end
    
    w = (-1).^(0:n-1)'.*[0.5;ones(n-2,1);0.5]; % Vector is 1/2 -1 1 ... -1 1 -1/2 (this is 1./c in Tref)
    
    j = (1:n+1:n^2);    % DNS: Vector between 1 and n^2 split into n-1 equal steps, i.e. of length n.
	k = ones(n,1);      % DNS: Vector of 1s length n
    Dw = w(:,k);        % DNS: Dw is a matrix size n by n, with w in every column, equivalent to repmat(w,1,10)
    Dw = Dw' ./ Dw - eye(n);    % DNS: Without the -eye(n), this is the same as (c*(1./c)') in Tref
    Dx = x(:,k);        % DNS: Dx is a matrix size n by n, with x in every column, equivalent to repmat(x,1,10)
    Dx = Dx - Dx' + eye(n);     % DNS: This gets the same Dx as dX+(eye(N)) in Tref
    D1 = Dw./Dx;        % DNS: The off diagonals now agree with Tref for D1 here, but with zeros on the diagonal, Tref has 1s
    D1(j) = 0;          % DNS: Make sure diagonals are zero (they should be already)
    d1 = -D1*k;         % DNS: Same as adding up the rows of D1 (with a -ve sign out the front), same as -sum(D1,2)
    D1(j) = d1;         % DNS: Put the d1 on the diagonals. So then this and Tref agree.
    Diff.Dx = D1;
    
    if (order == 1), return; end
    D2 = 2*(Dw.*d1(:,k)-D1)./Dx;      % DNS: 2 is the number of diffs (ell in Weid/Reddy). Rest is comparable to their code
    D2(j) = 0;                        % DNS: Make diagonals of D2 zero
    d2 = -D2*k;                       % DNS: Correct the diagonals, same as -sum(D2,2), or indeed -sum(D2')
    D2(j) = d2;                       % DNS: Put them back into D2
    
    Diff.DDx = D2;
    
    if (order == 2), return; end
    
    D3 = 3*(Dw.*d2(:,k)-D2)./Dx;      % DNS: D3 and D4 similarly
    D3(j) = 0; d3 = -D3*k; D3(j) = d3;
    
    Diff.DDDx = D3;
    
    if (order == 3), return; end
    
    D4 = 4*(Dw.*d3(:,k)-D3)./Dx;
    D4(j) = 0; d4 = -D4*k; D4(j) = d4;
    
    Diff.DDDDx = D4;
    
end