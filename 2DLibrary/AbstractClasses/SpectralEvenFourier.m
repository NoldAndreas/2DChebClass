classdef (Abstract) SpectralEvenFourier < SpectralFourier
        
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = SpectralEvenFourier(N1,N2)
             this@SpectralFourier(N1,2*N2-2);                                    
             
             this.N2     = N2;
             this.Pts.x2 = this.Pts.x2(1:N2);  % = [0,1/(2*N2-1),2/(2*N2-1),...,0.5]
        end      
    end
    
    %**********************************************
    %************   Initializations    *************
    %**********************************************
    methods (Access = public) 
        
         function int = ComputeIntegrationVector(this)
            [h1,dy1]        = this.PhysSpace1(this.Pts.x1);
            [h1,dy2]        = this.PhysSpace2(this.Pts.x2);             
            
            [h1,wInt1] = ClenCurtFlip(this.N1-1); 
            wInt2     = 0.5*[0.5,ones(1,this.N2-2),0.5]/(this.N2-1);
            %1/N2 weight of integration in second direction
            int       = kron(dy1'.*wInt1,(dy2').*wInt2);            
            this.Int  = int;             
            
         end         
         function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)
            if(nargin < 5)
                saveBool = false;
            end
            if(nargin < 4)
                fullInterpGrid = false;
            end                
            %note that the Interpolation function in SpectralFourier uses 
            %the FFTMatrix of SpectralEvenFourier, hence this should be
            %fine.
            Interp = ComputeInterpolationMatrix@SpectralFourier(this,interp1,interp2,fullInterpGrid,saveBool,2*this.N2-2);                 
        end
        function [fftMatrix,IfftMatrix] = getFFTMatrices(this,N1,N2)

            n = ((1:N2) -1);
            n2 = N2/2+1;

            fftMatrix = exp( -1i * 2*pi *(n'*n) /(N2));
            
            IM = [zeros(n2-2,1),fliplr(eye(n2-2)),zeros(n2-2,1)];

            fftMatrix = fftMatrix*[(eye(n2));IM];

            IfftMatrix = 0;%1/N2 * exp( 1i * 2*pi *(n'*n) /N2);

            fftMatrix = kron(eye(N1),fftMatrix);
            IfftMatrix = 0;%kron(eye(N1),IfftMatrix);

        end        
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
           
            [y1_kv,dy1] = PhysSpace1(this,x1);
            [y2_kv,dy2] = PhysSpace2(this,x2);
            
            if(nargout > 2)
                n        = length(x1);
                J        = zeros(n,2,2);
                J(:,1,1) = dy1;
                J(:,2,2) = dy2;                
            end
            
            if(nargout > 3)
                dH1 = 0;
                dH2 = 0;
                disp('SpectralFourierClass:Comp_to_Phys: dH1/2 NOT YET IMPLEMENTED');
            end
            
        end
        function [x1,x2] = CompSpace(this,y1,y2)
            x1 = CompSpace1(this,y1);
            x2 = CompSpace2(this,y2);
        end
    end
    
    %**********************************************
    %************   Mapping functions *************
    %**********************************************
    methods (Abstract = true,Access = public)         
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(x);
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(x); %from [0 0.5]->Physical Space
         x = CompSpace1(y);
         x = CompSpace2(y);
    end  
    
end


