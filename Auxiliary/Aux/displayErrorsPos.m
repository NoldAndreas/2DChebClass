
function err = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,cart)
    err = struct();

    [err.interp,err.interp_i] = max(abs(VP - vplot));
        
    str = ['Error of Interpolation:',num2str(err.interp),'\n'];
    if(err.interp > 10^-1)
        cprintf('red',str);
    elseif(err.interp <= 10^(-5))
        cprintf(str);
    else
        cprintf('magenta',str);
    end                       	
    
    if(nargin < 5)
        return;
    end
    
    if(isfield(Diff,'Dy1') && isfield(Vdiff,'dy1'))
        [err.dy1,err.dy1_i] = GetMaxAndPos(Diff.Dy1,Vdiff.dy1,'d/dy1');        
    end
    if(isfield(Diff,'Dy2') && isfield(Vdiff,'dy2'))
        [err.dy2,err.dy2_i] = GetMaxAndPos(Diff.Dy2,Vdiff.dy2,'d/dy2');                
    end
    if(isfield(Diff,'DDy1') && isfield(Vdiff,'ddy1'))
        [err.ddy1,err.ddy1_i] = GetMaxAndPos(Diff.DDy1,Vdiff.ddy1,'d^2/dy1^2');                        
    end
    if(isfield(Diff,'DDy2') && isfield(Vdiff,'ddy2'))
        [err.ddy2,err.ddy2_i] = GetMaxAndPos(Diff.DDy2,Vdiff.ddy2,'d^2/dy2^2');                        
    end
    if(isfield(Diff,'Dy1Dy2') && isfield(Vdiff,'dy1dy2'))
        [err.dy1dy2,err.dy1dy2_i] = GetMaxAndPos(Diff.Dy1Dy2,Vdiff.dy1dy2,'d^2/dy1dy2');                        
    end
    if(isfield(Diff,'Lap') && isfield(Vdiff,'Lap'))
        [err.Lap,err.Lap_i] = GetMaxAndPos(Diff.Lap,Vdiff.Lap,'Laplace');                        
    end
    if(isfield(Diff,'DDDy1') && isfield(Vdiff,'dddy1'))
        [err.dddy1,err.dddy1_i]= max(abs(Diff.DDDy1*V - Vdiff.dddy1));
        display(['Error of d^3/dy1^3:',num2str(err.dddy1),' , max at y1 = ',num2str(Pts.y1_kv(err.dddy1_i))]);        
    end
    if(isfield(Diff,'DDDDy1') && isfield(Vdiff,'ddddy1'))
        [err.ddddy1,err.ddddy1_i] = max(abs(Diff.DDDDy1*V - Vdiff.ddddy1));
        display(['Error of d^4/dy1^4:',num2str(err.ddddy1),' , max at y1 = ',num2str(Pts.y1_kv(err.ddddy1_i))]);        
    end

    if((nargin == 6) && (strcmp(cart,'cart')))
        if(isfield(Diff,'Lap') && isfield(Vdiff,'ddy1') && isfield(Vdiff,'ddy2'))
            h = Diff.Lap*V - (Vdiff.ddy1 + Vdiff.ddy2);
            [err.Lap,err.Lap_i] = max(abs(h));
            display(['Error of Lap:',num2str(err.Lap),' , max at y1 = ',num2str(Pts.y1_kv(err.Lap_i))]);
        end        
        if(isfield(Diff,'grad') && isfield(Vdiff,'dy1') && isfield(Vdiff,'dy2'))
            h = Diff.grad*V - [Vdiff.dy1 ; Vdiff.dy2];
            [err.Grad,err.Grad_i] = max(abs(h));
            display(['Error of grad:',num2str(err.Grad),' , max at y1 = ',num2str(Pts.y1_kv(err.Grad_i))]);
        end        
    end
    
    function [error,ind] = GetMaxAndPos(dOp,dCheck,text)
        [error,ind] = max(abs(dOp*V - dCheck));
        
        str = ['Error of ',text,':',num2str(error),...
               ' , max at y1 = ',num2str(Pts.y1_kv(ind)),...
               ' , max at y2 = ',num2str(Pts.y2_kv(ind)),'\n'];
        if(error > 10^-1)
            cprintf('red',str);
        elseif(error <= 10^(-5))
            cprintf(str);
        else
            cprintf('magenta',str);
        end                
    end
end