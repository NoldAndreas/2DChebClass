function I  = GetIndicesTrapezium(this)
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
% - indices of 'finite', 'infinite' boundaries
% - normal vectors: normalFinite, normalInfinite
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
        
        finite1 = ( (leftFinite & left) | (rightFinite & right) );
        finite2 = ( (topFinite & top) | (bottomFinite & bottom) );        
        finite  = (finite1 | finite2);
               
        infinite = ( (~leftFinite & left) | (~rightFinite & right) ...
                   | (~topFinite & top) | (~bottomFinite & bottom) );
                              
        corners = ((right & top) | (top & left) | (left & bottom) | (bottom & right));
        
        Z             = zeros(length(x1_kv));

        nx1Left       = Z;
        nx1Right      = Z;
        nx2Left       = Z;
        nx2Right      = Z;
        
        nx2Top        = Z;
        nx2Bottom     = Z;
         
%         n1 = zeros(N2,1);
%         n2 = zeros(N1,1);
%         
%         for i=1:N2 
%         %n1(i) = ?;
%         end
%         for i=1:N1
%         %n2(i) = ?;
%         end
%         %n1 = n1./sqrt(n1.^2);
%         %n2 = n2./sqrt(n2.^2);
        
        n1Temp = 1;
        n2Temp = (this.L2 - this.L1)/(2*this.h);
        normN = sqrt(n1Temp^2+n2Temp^2);
        
        n1 = n1Temp/normN;
        n2 = n2Temp/normN;

        nx1Left(left,left)     = -n1*speye(sum(left));
        nx2Left(left,left)     = -n2*speye(sum(left));
        nx1Right(right,right)  = n1*speye(sum(right));
        nx2Right(right,right)  = -n2*speye(sum(right));
        
        nx2Top(top,top)        = speye(sum(top));
        nx2Bottom(bottom,bottom)  = -speye(sum(bottom));
        
        nf1 = leftFinite*nx1Left(finite,:) + rightFinite*nx1Right(finite,:);
        nf2 = leftFinite*nx2Left(finite,:) + rightFinite*nx2Right(finite,:) + ...
              topFinite*nx2Top(finite,:)   + bottomFinite*nx2Bottom(finite,:);
        
        normalFinite = sparse([nf1 nf2] );
                    
        normalFinite1 = sparse( nf1 );
        normalFinite2 = sparse( nf2 );
                    
        normalInfinite = sparse( [((~leftFinite)*nx1Left(infinite,:) + (~rightFinite)*nx1Right(infinite,:)) ...
                        ((~topFinite)*nx2Top(infinite,:) + (~bottomFinite)*nx2Bottom(infinite,:))] );
        
        I = struct('left',left,'right',right,'bottom',bottom,'top',top,...
            'bound',bound,...
            'normalLeft',  sparse( [nx1Left(left,:) nx2Left(left,:) ] ),...
            'normalRight', sparse( [nx1Right(right,:) nx2Right(right,:)] ),...  % this had a nx2Left instead of nx2Right
            'normalLeft1',  sparse( nx1Left(left,:) ),...
            'normalRight1', sparse( nx1Right(right,:) ),...
            'normalLeft2',  sparse( nx2Left(left,:) ),...
            'normalRight2', sparse( nx2Right(right,:) ),...
            'normalTop',   sparse( [Z(top,:)  nx2Top(top,:)] ),...
            'normalBottom',sparse( [Z(bottom,:)  nx2Bottom(bottom,:)] ),...
            'normal',sparse( [(nx1Left(bound,:)+nx1Right(bound,:)) (nx2Left(bound,:) + nx2Right(bound,:) + nx2Top(bound,:) + nx2Bottom(bound,:))] ),...
            'corners',corners, ...
            'finite',finite,'infinite',infinite, ...
            'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
            'finite1',finite1,'finite2',finite2,...
            'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);
        
    end