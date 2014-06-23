classdef (Abstract) M1SpectralSpectral < Shape
    
    properties
        borderRight,borderTop,borderLeft,borderBottom
    end
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = M1SpectralSpectral(N1,N2)
             this@Shape(N1,N2);
             
             this.Pts.x1 = ClenCurtFlip(N1-1); 
             this.Pts.x2 = ClenCurtFlip(N2-1);           
             this.M      = N1*N2;
        end      
    end
    
    %**********************************************
    %************   Initializations    *************
    %**********************************************
    methods (Access = public) 

        function Int    = ComputeIntegrationVector(this,t1Odd,t2Odd)
            if((nargin >= 2) && t1Odd)
                wInt1O  = pi*[0.5 ones(1,this.N1-2) 0.5]/(this.N1-1);%sin(theta)!! 
                t       = pi*(0:(this.N1-1))/(this.N1-1);
                wInt1   = wInt1O.*sin(t);
            else
                [h_1,wInt1]  = ClenCurtFlip(this.N1-1);  
            end
            if((nargin >= 3) && t2Odd)
                wInt2      = [0.5 ones(1,this.N2-2) 0.5]/(this.N2-2);
            else
                [h_1,wInt2]  = ClenCurtFlip(this.N2-1);  
            end             
            
            [h_1,h_2,J]    = PhysSpace(this,this.Pts.x1_kv,this.Pts.x2_kv);            
            det        = DetMatr(J);
            det        = abs(det);
            det(det == inf) = 0; %Here, we assumed that integration is well-defined -> inf
            det(isnan(det)) = 0; %Here, we assumed that integration is well-defined -> inf
            Int        = kron(wInt1,wInt2).*(det');
            this.Int = Int;
        end        
        function Diff   = ComputeDifferentiationMatrix(this)
            
            Diff1 = barychebdiff(this.Pts.x1,2);
            Diff2 = barychebdiff(this.Pts.x2,2);    
            [h1,h2,J,dH1,dH2]    = PhysSpace(this,this.Pts.x1_kv,this.Pts.x2_kv);                        
            Diff      = PhysicalDerivativeJH(this.Pts,Diff1,Diff2,J,dH1,dH2);
            this.Diff = Diff;
        end              
        function [Interp,Interp1,Interp2] = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)            
            % Spectral Interpolation
            % Interpolates from the usual Chebyshev grid to the plotting grid given
            % above..    
            [Interp1,Interp2] = ComputeInterpolationMatrix12(this,interp1,interp2);
            %Interp1 = barychebevalMatrix(this.Pts.x1,interp1);  
            %Interp2 = barychebevalMatrix(this.Pts.x2,interp2);

            % kron form of the two interpolations                                                            
            Interp = struct('InterPol',kron(Interp1,Interp2),...
                            'Nplot1',length(interp1),...
                            'Nplot2',length(interp2),...
                            'N1',this.N1,'N2',this.N2);
                        
            if((nargin >= 4) && fullInterpGrid)
                x1 = kron(interp1,ones(size(interp2)));
                x2 = kron(ones(size(interp1)),interp2);
                [Interp.pts1,Interp.pts2] = PhysSpace(this,x1,x2);
            end
            if((nargin >= 5) && saveBool)
                this.Interp = Interp;
            end                        
        end
        function [I1,I2] = ComputeInterpolationMatrix12(this,interp1,interp2)
            I1 = barychebevalMatrix(this.Pts.x1,interp1);  
            I2 = barychebevalMatrix(this.Pts.x2,interp2);
        end
        function Ind    = ComputeIndices(this)
            Ind      = GetIndicesBox(this);
            this.Ind = Ind;
        end                 
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            N1  = this.N1;  N2  = this.N2;
            Pts = this.Pts;
            Int = this.Int;  % 1 x N1*N2

            if(useDistance)
                fPTemp = f(GetDistance(this,Pts.y1_kv,Pts.y2_kv));
            else
                fPTemp = f(Pts.y1_kv,Pts.y2_kv);
            end
            
            fDim = size(fPTemp);
            
            nElts = prod(fDim(2:end));
            
            IntT = Int.';  % N1*N2 x 1
            
            IntT = IntT(:,ones(1,nElts)); % N1*N2 x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([N1*N2,N1*N2,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:(N1*N2) 
                if(useDistance)
                    fP          = f(GetDistance(this,Pts.y1_kv(i) - Pts.y1_kv,Pts.y2_kv(i) - Pts.y2_kv));
                else
                    fP          = f(Pts.y1_kv(i) - Pts.y1_kv,Pts.y2_kv(i) - Pts.y2_kv);
                end
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;
            
            if((nargin >= 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end                
        function SetUpBorders(this,N)
            
            if(length(N) == 1)
                N = ones(4,1)*N;
            end
            
            this.borderRight = Border(N(1),this,'right',this.polar);
            this.borderRight.ComputeIntegrationVector();
            
            this.borderTop = Border(N(2),this,'top',this.polar);
            this.borderTop.ComputeIntegrationVector();
            
            this.borderLeft = Border(N(3),this,'left',this.polar);
            this.borderLeft.ComputeIntegrationVector();
            
            this.borderBottom = Border(N(4),this,'bottom',this.polar);
            this.borderBottom.ComputeIntegrationVector();
            
        end
        function int = IntFluxThroughDomain(this,N)
            
           SetUpBorders(this,N);
           int = this.borderRight.IntNormal + this.borderTop.IntNormal + ...
                 this.borderLeft.IntNormal + this.borderBottom.IntNormal;
        end        
        function PlotBorders(this)
            this.borderRight.PlotPath(this.polar);
            this.borderTop.PlotPath(this.polar);
            this.borderLeft.PlotPath(this.polar);
            this.borderBottom.PlotPath(this.polar);
        end
        
    end
    
end


