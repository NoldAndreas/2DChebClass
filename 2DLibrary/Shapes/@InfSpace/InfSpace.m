classdef InfSpace < SpectralSpectral

    properties                
        AD = [];
        L1,L2
        y10 = 0;
        y20 = 0;
    end
    
    methods        
        function this = InfSpace(Geometry)
            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));
                   
            this.polar = 'cart';
            this.L1    = Geometry.L1;
            this.L2    = Geometry.L2;
            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end
            if(isfield(Geometry,'y20'))
                this.y20 = Geometry.y20;
            end            
            InitializationPts(this);                        
        end                         
    end
    
    methods (Access = public)
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)                        
            [y1,dy1,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);
            y1 = y1 + this.y10;
        end
        function x1 = CompSpace1(this,y1)            
            x1  = InvSqrtMap(y1-this.y10,this.L1,inf);
        end
        function [y2,dy2,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                        
            [y2,dy2,dx,ddx,dddx,ddddx] = SqrtMap(x2,this.L2,inf);
            y2 = y2 + this.y20;
        end
        function x2 = CompSpace2(this,y2)            
            x2  = InvSqrtMap(y2-this.y20,this.L2,inf);
        end
         
        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams,parent)
            %{
            Strategy:
            We wish to compute \int_D(x) f(x-y) g(y) dy, for some domain D
            (e.g. disc, annulus) centred at x.
            Let z = y-x then dy = dz, y = x + z and 
            y \in D(x) => z + x \in D(x) or z \in D(0) and we get
            \int_D(0) f(-z) g(x+z) dz
            
            In matlab this is
            M_conv = Int_D(0) .* f(-z) .* IP_D(x)
            %}
            
            if((nargin>=4) && parent)
                M_conv = ComputeConvolutionMatrix@M1SpectralSpectral(this,f,[]);
                return;
            end
            
            
            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end
            
            disp('Computing Convolution matrices...'); 

            n1 = this.N1;   n2 = this.N2;

            if(isfield(shapeParams,'R'))
                disp('Disc');
                subShape = Disc(shapeParams);
                rRange = (-1:0.05:1)';   
            elseif(isfield(shapeParams,'RMin'))
                if(isfield(shapeParams,'RMax'))
                    disp('Annulus');
                    subShape = Annulus(shapeParams);
                    rRange = (-1:0.05:1)';
                else
                    disp('InfAnnulus');
                    subShape = InfAnnulus(shapeParams);
                    rRange = (-1:0.05:0.9)';
                end
            else
                disp('InfDisc');
                subShape = InfDisc(shapeParams);
                rRange = (-0.9:0.05:0.9)';
            end
            
            % Integration and interpolation on D(0)
            Int_i    = subShape.ComputeIntegrationVector();
            subShape.ComputeInterpolationMatrix(rRange,(0:0.02:1)',true,true);
            
            % f(-z)
            subShapeCartPts = Pol2CartPts(subShape.Pts);
            invSubShapePts = invertPts(subShapeCartPts,'cart');
            
            if(useDistance)
                fP       = f(this.GetDistance(invSubShapePts.y1_kv,invSubShapePts.y2_kv));
            else
                fP       = f(invSubShapePts.y1_kv,invSubShapePts.y2_kv);
            end
 
            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end
            
            nf = size(fP,2);
            
            IntRep = Int_i(ones(1,nf),:);
            
            M_conv = zeros(n1*n2,n1*n2,nf);  
                        
            figure
            subShape.doPlots(fP);
            title('Interpolation of Convolution function'); drawnow;

            
            y1_kv=this.Pts.y1_kv;
            y2_kv=this.Pts.y2_kv;
            
            hw = waitbar(0,'Computing Convolution matrices...');

            for i=1:(n1*n2)
                
                waitbar(i/(n1*n2),hw);
                
                Y10 = y1_kv(i); 
                Y20 = y2_kv(i);

                % Compute points in D(x) = (Y10,Y20))
                subConvShapePts = shiftPointsCart(subShapeCartPts,Y10,Y20);
                % and the interpolation onto D(x) from the full domain
                IP          = SubShapePts(this,subConvShapePts);
                
                IntFTemp = (IntRep.*fP.')*IP;
                M_conv(i,:,:) = IntFTemp.';
                
            end
            
            close(hw);
            
            M_conv(isnan(M_conv)) = 0;

            if(nDim>2)
                nfTemp = size(fP_old);
                nfTemp = nfTemp(2:end);
                M_conv = reshape(M_conv,[ n1*n2, n1*n2, nfTemp ]);
            end
            
        end
        
        function test(this)            
            %this.ComputeConvolutionMatrix_Test;            
            %this.ComputeConvolutionMatrix_Test_Short_1; % works for pointwise
            this.ComputeConvolutionMatrix_Test_Short_2;  % works for both
        end
        
    end 
    
    methods (Access = private)
        
        ComputeConvolutionMatrix_Test(this);
        testAll(this);            
        
    end
end