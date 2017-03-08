
classdef Circle < FourierPath
    properties 
        R
        Origin = [0,0]
    end
    
    methods
        function this = Circle(Geometry)
            this@FourierPath(Geometry.N,'polar');            
            this.R      = Geometry.R;       
            this.polar  = 'polar';
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            InitializationPts(this);
        end
        
        %***************************************************************
        % Mapping functions:
        %***************************************************************
        function [y1 ,y2 , dy1_dt , dy2_dt ]  = f_path(this,t)        
            y1      = this.R*ones(size(t));   y2      = 2*pi*t;
            dy1_dt  = zeros(size(t));         dy2_dt  = 2*pi*ones(size(t));
        end       
        function [int,length] = ComputeIntegrationVector(this)
            int    = ComputeIntegrationVector@FourierPath(this);
            length = 2*pi*this.R;
            disp(['Circle: Error of integration of length (ratio): ',...
                                    num2str(1-sum(this.IntSc)/length)]);
        end
        function [FFTMatrix,IFFTMatrix] = getFFTMatrices(this)

            N=this.N;
            
            n = (1:N) -1;
            k = (1:N) -1;

            FFTMatrix = exp( -1i * 2*pi *n'*k /N);
            IFFTMatrix = 1/N * exp( 1i * 2*pi *n'*k /N);
        end              
        function Interp = ComputeInterpolationMatrix(this,interp)
            
            N = this.N;
            
            LL    = 2*pi;
            M     = floor((N+1)/2);

            %%FFT Interpolation
            n21   = 1:(M-1);       % first half of points, excluding zero
            n22   = (M+1):(N-1);  % second half

            % This uses that NT+1 is odd (ie NT is even)
            phi      = interp*2*pi;    
            InterpTemp  = [exp(2*pi*1i*phi/LL*0) , exp(2*pi*1i*phi/LL*n21) , cos(2*pi*M*phi/LL) , exp(2*pi*1i*phi/LL*(n22-N))]/N;
        
            interp_y = interp;
            
            [FFTMatrix,h1s] = getFFTMatrices(this);
            
            InterPol = real(InterpTemp*FFTMatrix);
            
            Interp = struct('InterPol',InterPol,...
                'interp',interp_y,...             
                'Nplot',length(interp),...
                'N',this.N);                   
            
        end
        function PlotFunction(this,f,NPlot)

            if(nargin<3)
                NPlot=100;
            end
            
            t = this.Pts.y2_kv;
            r = ones(size(t))*this.R;
            [X,Y] = pol2cart(t,r);
            
            xPlot = (0:1/NPlot:1)';
            plotT = xPlot*2*pi;
            plotR = ones(size(plotT))*this.R;
            [plotX,plotY] = pol2cart(plotT,plotR);
            
            if(nargin<2)
                plot(plotX,plotY,'b');
                hold on
                plot(X,Y,'go');
            else            
                Interp = ComputeInterpolationMatrix(this,xPlot);
                F = f(t);
                plotF = Interp.InterPol*F;

                plot3(X,Y,F,'go');
                hold on;
                plot3(plotX,plotY,plotF,'b');            
            end  
        end
        

        function ptsCart = GetCartPts(this,pts_y1,pts_y2)
            
            if(nargin < 3)
                ptsCart  = GetCartPts@FourierPath(this);
            else
                ptsCart  = GetCartPts@FourierPath(this,pts_y1,pts_y2);                
            end
                        
            ptsCart.y1_kv = ptsCart.y1_kv + this.Origin(1);
            ptsCart.y2_kv = ptsCart.y2_kv + this.Origin(2);
        end
        
    end
   
end