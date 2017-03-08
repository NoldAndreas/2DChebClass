function acc = testIntegration(this,toCheck)

    
    if(nargin == 1)
        toCheck = '(y-y0).^2';
    end
    R    = this.R;
    int  = this.Int;
    
    th    = this.Pts.y1_kv;
    phi   = this.Pts.y2_kv;
    y     = R.*cos(th);

	y0     = 3;
    k      = 2*pi;
    theta1 = this.theta1;
    theta2 = this.theta2;
    acc    = 0;
    
    if(strcmp(toCheck,'area') || strcmp(toCheck,'all'))            
        analytical  = 2*this.R^2*pi*(cos(theta1)-cos(theta2));  
        numerical   = sum(int);
        acc         = abs(analytical-numerical);
        if(nargout==0)
            disp(['Error (absolute) for area:',num2str(acc)]);        
        end
    end
    if(strcmp(toCheck,'(y-y0).^2') || strcmp(toCheck,'all'))                                    
        analytical  = (2/3)*R^4*cos(theta1)^3*pi-2*R^3*cos(theta1)^2*y0*pi+2*R^2*cos(theta1)*y0^2*pi-(2/3)*R^4*cos(theta2)^3*pi+2*R^3*cos(theta2)^2*y0*pi-2*R^2*cos(theta2)*y0^2*pi;
        f           = (y-y0).^2;
        numerical   = sum(int.*f');
        acc         = acc + abs(analytical-numerical);
        if(nargout==0)
            disp(['Error (absolute) for (y-y0).^2:',num2str(abs(analytical-numerical))]);        
        end
    end
    if(strcmp(toCheck,'sin(y-y0)') || strcmp(toCheck,'all'))                                    
        analytical  = -2*R*(cos(k*(R*cos(theta1)-y0))-cos(k*(R*cos(theta2)-y0)))*pi/k;
        f           = sin(k*(y-y0));
        numerical   = sum(int.*f');
        acc         = acc + abs(analytical-numerical);
        if(nargout==0)
            disp(['Error (absolute) for sin(k*(y-y0).^2:',num2str(abs(analytical-numerical))]);        
        end
    end
    if(strcmp(toCheck,'exp(k*(y-y0))') || strcmp(toCheck,'all'))  
        k  = -1;
        y0 = 2;
        
        analytical  = 2*R*(exp(k*(R*cos(theta1)-y0))-exp(k*(R*cos(theta2)-y0)))*pi/k;
        f           = exp(k*(y-y0));
        numerical   = sum(int.*f');
        acc         = acc + abs(analytical-numerical);
        if(nargout==0)
            disp(['Error (absolute) for exp(k*(y-y0)):',num2str(abs(analytical-numerical))]);        
            disp(['Error (realative) for exp(k*(y-y0)):',num2str(abs(1-numerical/analytical))]);        
        end
    end    

end