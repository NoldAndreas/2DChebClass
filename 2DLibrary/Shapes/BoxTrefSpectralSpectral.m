classdef BoxTrefSpectralSpectral < TrefSpectralSpectral
    
    methods (Access = public) 
        function this = BoxTrefSpectralSpectral(N1,N2)
            this@TrefSpectralSpectral(N1,N2);            
            this.polar = 'cart';
            InitializationPts(this);
        end
                
        %********************************************************
        %****** Maps from [-1,1]x[-1,1] -> Physical Space
        % The global map will be
        % PhysSpace2Tref(Tref(x2,d,ep))
        % PhysSpace1Tref(x1)
        %********************************************************
        
        function [y,dy,ddy] = PhysSpace1Tref(this,x)
            y    = x;
            dy   = ones(size(x));
            ddy  = zeros(size(x));            
        end
        function [y,dy,ddy] = PhysSpace2Tref(this,x)
            y    = x;
            dy   = ones(size(x));
            ddy  = zeros(size(x));
        end
            
        function x = CompSpace1Tref(this,y)
            x = y;
        end
        function x = CompSpace2Tref(this,y)
            x = y;
        end
        
        %********************************************************
        
    end
    
end