
    function err = displayErrors(vplot,VP,V,Vdiff,Diff,cart)
        err = struct();
    
        display(['Error of Interpolation:',num2str(max(abs(vplot - VP)))]);
        if(isfield(Diff,'Dy1') && isfield(Vdiff,'dy1'))
            err.dy1 = max(abs(Diff.Dy1*V - Vdiff.dy1));
            dispError('Error of d/dy1:',err.dy1);
            %display(['Error of d/dy1:',num2str(err.dy1)]);
        end
        if(isfield(Diff,'Dy2') && isfield(Vdiff,'dy2'))
            err.dy2 = max(abs(Diff.Dy2*V - Vdiff.dy2));
            dispError('Error of d/dy2:',err.dy2);
            %display(['Error of d/dy2:',num2str(err.dy2)]);
        end
        if(isfield(Diff,'DDy1') && isfield(Vdiff,'ddy1'))
            err.ddy1 = max(abs(Diff.DDy1*V - Vdiff.ddy1));
            dispError('Error of d^2/dy1^2:',err.ddy1);
            %display(['Error of d^2/dy1^2:',num2str(err.ddy1)]);        
        end
        if(isfield(Diff,'DDy2') && isfield(Vdiff,'ddy2'))
            err.ddy2 = max(abs(Diff.DDy2*V - Vdiff.ddy2));
            dispError('Error of d^2/dy2^2:',err.ddy2);
            %display(['Error of d^2/dy2^2:',num2str(err.ddy2)]);
        end
        if(isfield(Diff,'Dy1Dy2') && isfield(Vdiff,'dy1dy2'))
            err.dy1dy2 = max(abs(Diff.Dy1Dy2*V - Vdiff.dy1dy2));
            dispError('Error of d^2/(dy1dy2):',err.dy1dy2);
            %display(['Error of d^2/(dy1dy2):',num2str(err.dy1dy2)]);
        end
        if(isfield(Diff,'DDDy1') && isfield(Vdiff,'dddy1'))
            err.dddy1 = max(abs(Diff.DDDy1*V - Vdiff.dddy1));
            dispError('Error of d^3/(dy1^3):',err.dddy1);
            %display(['Error of d^3/dy1^3:',num2str(err.dddy1)]);        
        end
        if(isfield(Diff,'DDDDy1') && isfield(Vdiff,'ddddy1'))
            err.ddddy1 = max(abs(Diff.DDDDy1*V - Vdiff.ddddy1));
            dispError('Error of d^2/(dy1^4):',err.ddddy1);
            %display(['Error of d^4/dy1^4:',num2str(err.ddddy1)]);        
        end
        
        if((nargin == 6) && (strcmp(cart,'cart')))
            if(isfield(Diff,'Lap') && isfield(Vdiff,'ddy1') && isfield(Vdiff,'ddy2'))
                h = Diff.Lap*V - (Vdiff.ddy1 + Vdiff.ddy2);
                err.Lap = max(abs(h));
                dispError('Error of Lap:',err.Lap);                
            end        
            if(isfield(Diff,'grad') && isfield(Vdiff,'dy1') && isfield(Vdiff,'dy2'))
                h = Diff.grad*V - [Vdiff.dy1 ; Vdiff.dy2];
                err.Grad = max(abs(h));
                dispError('Error of grad:',err.Grad);                
            end        
        end
        
        function dispError(str,val)
            if(val > 10^-3)
                cprintf('red',[str,num2str(val),'\n']);
            elseif(val <= 10^(-8))
                cprintf([str,num2str(val),'\n']);
            else
                cprintf('magenta',[str,num2str(val),'\n']);                
            end
        end
    end