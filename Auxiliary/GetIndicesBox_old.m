function I  = GetIndicesBox(x1,x2)
%************************************************************************
%I  = GetIndicesBox(x1,x2)
%INPUT:
% x1 - grid points in computational space in 1st variable, x1 in [-1 1]
% x2 - grid points in computational space in 2nd variable, x2 in [-1 1]
%OUTPUT: structure with ...
% - indices of 'left', 'right', 'bottom', 'top' boundaries
% - bound = includes left, right, top and bottom boundaries
% - normal vectors: normalLeft,normalRight,normalTop, normalBottom
% - normal = normal vectors for full boundary
% - corners 
%************************************************************************


        %Sparse Matrices are to be used
        N1 = length(x1);
        N2 = length(x2);

        x1_kv   = kron(x1,ones(size(x2)));
        x2_kv   = kron(ones(size(x1)),x2);
                
    
        left    = (x1_kv == min(x1_kv));
        right   = (x1_kv == max(x1_kv)); 
        bottom  = (x2_kv == min(x2_kv));
        top     = (x2_kv == max(x2_kv));
        
        bound   = (left | right | bottom | top);
        
        corners = ((right & top) | (top & left) | (left & bottom) | (bottom & right));
        
        Z             = zeros(N1*N2);
        
        nx1Left       = zeros(N1*N2);
        nx1Right      = zeros(N1*N2);
        
        nx2Top        = zeros(N1*N2);
        nx2Bottom     = zeros(N1*N2);
         
        nx1Left(left,left)     = -speye(sum(left));
        nx1Right(right,right)  = speye(sum(right));
        
        nx2Top(top,top)        = speye(sum(top));
        nx2Bottom(bottom,bottom)  = -speye(sum(bottom));
        
        I = struct('left',left,'right',right,'bottom',bottom,'top',top,...
            'bound',bound,...
            'normalLeft',  sparse( [nx1Left(left,:)  Z(left,:)] ),...
            'normalRight', sparse( [nx1Right(right,:)  Z(right,:)] ),...
            'normalTop',   sparse( [Z(top,:)  nx2Top(top,:)] ),...
            'normalBottom',sparse( [Z(bottom,:)  nx2Bottom(bottom,:)] ),...
            'normal',sparse( [(nx1Left(bound,:)+nx1Right(bound,:)) (nx2Top(bound,:) + nx2Bottom(bound,:))] ),...
            'corners',corners);
        
    end