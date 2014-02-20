classdef (Abstract) SpectralFourier < Shape
    
    properties (Access = public)
        FFTMatrix
        IFFTMatrix
    end
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = SpectralFourier(N1,N2)
             this@Shape(N1,N2);                                    
             if(mod(N2,2) == 1)
                 exception = MException('SpectralFourierClass:Constructor','N2 has to be even');
                 throw(exception);                 
             end                          
             
             this.Pts.x1 = ClenCurtFlip(N1-1); 
             this.Pts.x2 = (0:(N2-1))'/(N2);   
             
             this.Pts.N1 = N1;
             this.Pts.N2 = N2;
             
             [this.FFTMatrix,this.IFFTMatrix] = getFFTMatrices(this,N1,N2);     
        end      
    end
    
    %**********************************************
    %************   Initializations    *************
    %**********************************************
    methods (Access = public) 
        function this = InitializationPts(this)
            [y1,dy1]        = this.PhysSpace1(this.Pts.x1);
            [y2,dy2]        = this.PhysSpace2(this.Pts.x2); 
            this.Pts.y1_kv  = kron(y1,ones(size(y2)));
            this.Pts.y2_kv  = kron(ones(size(y1)),y2);
            this.M          = length(this.Pts.y1_kv);
        end
        function int = ComputeIntegrationVector(this)
            [h1,dy1]        = this.PhysSpace1(this.Pts.x1);
            [h1,dy2]        = this.PhysSpace2(this.Pts.x2);             
            
            [h1,wInt1] = ClenCurtFlip(this.N1-1); 
            %1/N2 weight of integration in second direction
            this.Int = kron(dy1'.*wInt1,(dy2').*ones(1,this.N2))/this.N2;            
            int = this.Int;
        end        
        function Diff = ComputeDifferentiationMatrix(this,pts)
            
            if(nargin == 1)
                pts = this.Pts;                
            end
            
            N2 = length(pts.x2);
            x1 = pts.x1;
            
            Diff1 = barychebdiff(x1,3);    
            Diff2 = struct();
            %[Dx1,DDx1,DDDx1]= barychebdiff(x1); 

            %Differentiation matrix in x2
            h         = 2*pi/N2;       
            column    = [0 .5*(-1).^(1:N2-1).*cot((1:N2-1)*h/2)];
            Diff2.Dx  = toeplitz(column,column([1 N2:-1:2]))*2*pi;% 2pi is needed to go back to the -1;1 - domain        


            h         = 2*pi/N2;       
            column2   = [(-(pi^2)/(3*h^2) - 1/6)  -.5*(-1).^(1:N2-1)./((sin((1:N2-1)*h/2)).^2)];
            Diff2.DDx = toeplitz(column2,column2([1 N2:-1:2]))*(2*pi)^2;        


            Sel = {'Dy1' ;'DDy1' ; 'Dy2'; 'DDy2' ;...
                  'Dy1Dy2'; 'Lap' ;'grad' ;'div';...
                  'gradDiv'; 'LapVec'}; 

            Diff = PhysicalDerivatives(this,pts,Sel,Diff1,Diff2,2);
            %Dx1,Dx2,DDx1,DDx2);        

        end              
        function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool,n2)
            Interp1     = barychebevalMatrix(this.Pts.x1,interp1); 

            %2nd variable: Fourier
            if(nargin == 6)
                N2 = n2;
            else
                N2 = this.N2;
            end
            LL    = 2*pi;
            M     = floor((N2+1)/2);

            %%FFT Interpolation
            n21   = 1:(M-1);       % first half of points, excluding zero
            n22   = (M+1):(N2-1);  % second half

            % This uses that NT+1 is odd (ie NT is even)
            phi      = interp2*2*pi;            
            Interp2  = [exp(2*pi*1i*phi/LL*0) , exp(2*pi*1i*phi/LL*n21) , cos(2*pi*M*phi/LL) , exp(2*pi*1i*phi/LL*(n22-N2))]/N2;

            InterPol = real(kron(Interp1,Interp2)*this.FFTMatrix);
            
            interp_y1 = this.PhysSpace1(interp1);
            interp_y2 = this.PhysSpace2(interp2);
            % kron form of the two interpolations                                                            
            Interp = struct('InterPol',InterPol,...
                            'interp1',interp_y1,'interp2',interp_y2,...             
                            'Nplot1',length(interp1),...
                            'Nplot2',length(interp2),...
                            'N1',this.N1,'N2',this.N2);                   
                        
            if((nargin >= 4) && fullInterpGrid)
                Interp.pts1 = kron(interp_y1,ones(size(interp_y2)));
                Interp.pts2 = kron(ones(size(interp_y1)),interp_y2);             
            end
            if((nargin >= 5) && saveBool)
                this.Interp = Interp;
            end                        
        end        
        function Ind = ComputeIndices(this)
            Ind      = GetIndicesBox(this.Pts.x1,this.Pts.x2);
            this.Ind = Ind;
        end    
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)

            N1 = this.N1;   N2 = this.N2;

            xs  = this.Pts.y1_kv; 
            ys  = this.Pts.y2_kv;            
            M_conv = zeros(N1*N2,N1*N2);

            for i=1:(N1*N2)
                M_conv(i,:)  = this.Int.*f(xs(i)-xs,ys(i)-ys)';
            end   
            if((nargin== 3) && saveBool)
                this.Conv = M_conv;
            end
            
        end
        
        function [FFTMatrix,IFFTMatrix] = getFFTMatrices(this,N1,N2)

            n = (1:N2) -1;
            k = (1:N2) -1;

            FFTMatrix = exp( -1i * 2*pi *n'*k /N2);
            IFFTMatrix = 1/N2 * exp( 1i * 2*pi *n'*k /N2);

            FFTMatrix = kron(eye(N1),FFTMatrix);
            IFFTMatrix = kron(eye(N1),IFFTMatrix);

        end

        %Maps
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
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(x);
         x = CompSpace1(y);
         x = CompSpace2(y);
    end  
    
end


