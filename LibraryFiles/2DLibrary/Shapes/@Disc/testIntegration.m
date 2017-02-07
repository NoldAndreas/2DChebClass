function acc = testIntegration(this,toCheck)
    
    if(nargin == 1)
        toCheck = '(x-x0).^2';
    end
    R    = this.R;
    int  = this.Int;
    
    r    = this.Pts.y1_kv;
    t    = this.Pts.y2_kv;
    x    = r.*cos(t);

	x0   = 3;
    
    if(this.sphere)                
        if(strcmp(toCheck,'area') || strcmp(toCheck,'all'))            
            analytical  = 4/3*pi*R^3;                
            numerical   = sum(int);
            acc         = abs(analytical-numerical);
            if(nargout==0)
                disp(['Error (absolute) for area:',num2str(acc)]);        
            end
        end
        if(strcmp(toCheck,'(x-x0).^2') || strcmp(toCheck,'all'))                                    
            analytical = pi*((4/3)*x0^2*R^3+(4/15)*R^5);
            f          = (x-x0).^2;
            numerical  = sum(int.*f');
            acc        = acc + abs(analytical-numerical);
            if(nargout==0)
                disp(['Error (absolute) for (x-x0).^2:',num2str(abs(analytical-numerical))]);        
            end
        end
    else
        if(strcmp(toCheck,'area') || strcmp(toCheck,'all'))                        
            analytical = pi*R^2;
            numerical  = sum(int);
            acc        = abs(analytical-numerical);
            if(nargout==0)
                disp(['Error (absolute) for area:',num2str(acc)]);
            end
        end
        if(strcmp(toCheck,'(x-x0).^2') || strcmp(toCheck,'all'))                                                
            analytical = R^2*x0^2*pi+(1/4)*R^4*pi;                                
            f          = (x-x0).^2;
            numerical  = sum(int.*f');
            acc        = acc+abs(analytical-numerical);
            if(nargout==0)
                disp(['Error (absolute) for (x-x0).^2:',num2str(abs(analytical-numerical))]);        
            end
        end
    end        

end