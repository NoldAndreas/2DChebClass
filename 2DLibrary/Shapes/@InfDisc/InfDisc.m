classdef InfDisc < Polar_SpectralFourier

    properties        
        L
    end
    
    methods        
        function this = InfDisc(Geometry)
            this@Polar_SpectralFourier(Geometry.N(1),Geometry.N(2));
            this.L = Geometry.L;            
            
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            
            InitializationPts(this);
        end                         
    end
    
    methods (Access = public)        
        
        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            [r,dr,dx,ddx,dddx,ddddx] = SqrtMap(x,this.L,inf);             
        end
        function xf = CompSpace1(this,r)
            xf = InvSqrtMap(r,this.L,inf);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,xT)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end
        function xf = CompSpace2(this,th)
            xf = th/(2*pi);
        end        
        
        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams,parent)
            if((nargin == 4) && parent)
                M_conv = ComputeConvolutionMatrix@Polar_SpectralFourier(this,f);
                return;
            end
            
            %{
            Strategy:
            We wish to compute \int_D(x) f(x-y) g(y) dy, for some domain D
            (e.g. disc, annulus) centred at x.
            Let z = y-x then dt = dx, y = x + z and 
            y \in D(x) => z + x \in D(x) or z \in D(0) and we get
            \int_D(0) f(-z) g(x+z) dz
            
            In matlab this is
            M_conv = Int_D(0) .* f(-z) .* IP_D(x)
            %}
            
            disp('Computing Convolution matrices...'); 

            n1 = this.N1;   n2 = this.N2;
            M_conv = zeros(n1*n2,n1*n2);  
            
            if(isfield(shapeParams,'R'))
                subShape = Disc(shapeParams);
                rRange = (-1:0.05:1)';
            elseif(isfield(shapeParams,'RMin'))
                if(isfield(shapeParams,'RMax'))
                    subShape = Annulus(shapeParams);
                    rRange = (-1:0.05:1)';
                else
                    subShape = InfAnnulus(shapeParams);
                    rRange = (-1:0.05:0.9)';
                end
            else
                subShape = InfDisc(shapeParams);
                rRange = (-0.9:0.05:0.9)';
            end
            
            % Integration and interpolation on D(0)
            Int_i    = subShape.ComputeIntegrationVector();
            Interp   = subShape.ComputeInterpolationMatrix(rRange,(0:0.02:1)',true,false);
            
            % f(-z)
            invSubShapePts = invertPts(subShape.Pts,'polar');
            fP       = f(GetDistance(this,invSubShapePts.y1_kv,invSubShapePts.y2_kv));

            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end
            
            nf = size(fP,2);
            
            IntRep = Int_i(ones(1,nf),:);
            
            M_conv = zeros(n1*n2,n1*n2,nf);  

           figure
           subShape.Interp=Interp;
           subShape.plot(fP,'SC');
%            title('Interpolation of Convolution function'); drawnow;

            y1_kv=this.Pts.y1_kv;
            y2_kv=this.Pts.y2_kv;
            
            hw = waitbar(0,'Computing Convolution matrices...');
            
            for i=1:(n1*n2)
                
                waitbar(i/(n1*n2),hw);
                
                y10 = y1_kv(i); 
                y20 = y2_kv(i);

                % Compute points in D(x) = (y10,y20))
                subConvShapePts = shiftPointsPolar(subShape.Pts,y10,y20,'polar');
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
        
        [M_conv,M_conv_New] = ComputeConvolutionMatrix_Test(this,shapeParams)
        
    end
end