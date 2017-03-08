function [y,Int,Diff,Diff2] = InitAnalysisGridY(this,yInt,m,CHEB)

   if((nargin > 3) && strcmp(CHEB,'CHEB'))
        [x,w]  = ClenCurtFlip(m);
        y      = yInt(1) + (x+1)/2*(yInt(2)-yInt(1));
        Int    = w*(yInt(2)-yInt(1))/2;
        D      = barychebdiff(y,2);
        Diff   = D.Dx;
        Diff2  = D.DDx;
    else
        y      = yInt(1) + (0:(m-1))'/(m-1)*(yInt(2)-yInt(1));        
        Int    = ([y(2:end)-y(1:end-1);0]/2 + [0;y(2:end)-y(1:end-1)])'/2;
                          
        e1     = [1,0.5*ones(1,m-2)];
        e2     = [0.5*ones(1,m-2),1];
        D      = diag(e1,1)-diag(e2,-1);
        D(1,1) = -1;  D(m,m) = 1;
        Diff   = D*(m-1)/(yInt(2)-yInt(1));          
    end
    
end