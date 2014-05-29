function I  = GetIndicesBox_New(this)
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

        x1_kv = this.Pts.x1_kv;
        x2_kv = this.Pts.x2_kv;
        
        y1_kv = this.Pts.y1_kv;
        y2_kv = this.Pts.y2_kv;
        
        left    = (x1_kv == min(x1_kv));
        right   = (x1_kv == max(x1_kv)); 
        bottom  = (x2_kv == min(x2_kv));
        top     = (x2_kv == max(x2_kv));
                
        bound   = (left | right | bottom | top);
        
        leftFinite   = any(isfinite(y1_kv(left)));
        rightFinite  = any(isfinite(y1_kv(right)));
        topFinite    = any(isfinite(y2_kv(top)));
        bottomFinite = any(isfinite(y2_kv(bottom)));
        
        finite = ( (leftFinite & left) | (rightFinite & right) ...
                   | (topFinite & top) | (bottomFinite & bottom) );
        
        infinite = ( (~leftFinite & left) | (~rightFinite & right) ...
                   | (~topFinite & top) | (~bottomFinite & bottom) );
               
               
        corners = ((right & top) | (top & left) | (left & bottom) | (bottom & right));
        
        Z             = zeros(length(x1_kv));

        nx1Left       = Z;
        nx1Right      = Z;
        
        nx2Top        = Z;
        nx2Bottom     = Z;
         
        nx1Left(left,left)     = -speye(sum(left));
        nx1Right(right,right)  = speye(sum(right));
        
        nx2Top(top,top)        = speye(sum(top));
        nx2Bottom(bottom,bottom)  = -speye(sum(bottom));
        
        normalFinite = sparse( [(leftFinite*nx1Left(bound,:) + rightFinite*nx1Right(bound,:)) ...
                        (topFinite*nx2Top(bound,:) + bottomFinite*nx2Bottom(bound,:))] );
        normalInfinite = sparse( [((~leftFinite)*nx1Left(bound,:) + (~rightFinite)*nx1Right(bound,:)) ...
                        ((~topFinite)*nx2Top(bound,:) + (~bottomFinite)*nx2Bottom(bound,:))] );
        
        I = struct('left',left,'right',right,'bottom',bottom,'top',top,...
            'bound',bound,...
            'normalLeft',  sparse( [nx1Left(left,:)  Z(left,:)] ),...
            'normalRight', sparse( [nx1Right(right,:)  Z(right,:)] ),...
            'normalTop',   sparse( [Z(top,:)  nx2Top(top,:)] ),...
            'normalBottom',sparse( [Z(bottom,:)  nx2Bottom(bottom,:)] ),...
            'normal',sparse( [(nx1Left(bound,:)+nx1Right(bound,:)) (nx2Top(bound,:) + nx2Bottom(bound,:))] ),...
            'corners',corners, ...
            'finite',finite,'infinite',infinite, ...
            'normalFinite',normalFinite,'normalInfinite',normalInfinite);
        
    end