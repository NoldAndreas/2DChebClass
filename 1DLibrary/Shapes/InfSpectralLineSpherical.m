classdef InfSpectralLineSpherical < Spectral
    
    properties
        L
        y0=0;
    end
    
    methods
        function this = InfSpectralLineSpherical(Geometry)
            this@Spectral(Geometry.N);
            
            this.L = Geometry.L;

            this.polar = 'polar';
            
            InitializationPts(this);  
        end
        
    end
       
   methods (Access = public)
        
       function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace(this,x)
            [y,dy,dx,ddx,dddx,ddddx] = SqrtMap(x,this.L,inf);
        end
        
        function xf = CompSpace(this,y)
            xf  = InvSqrtMap(y-this.y0,this.L,inf);
        end       
        
         function w = IntegrateRegion(this,weight,shapeParams)
            if(~isfield(shapeParams,'yMin'))
                shapeParams.yMin = 0;
            end
            if(~isfield(shapeParams,'yMax'))
                shapeParams.yMax = inf;
            end
            
            if(~isfield(shapeParams,'L'))
                shapeParams.L = this.L;
            end
            if(~isfield(shapeParams,'N'))
                shapeParams.N = this.N;
            end
            
            if(shapeParams.yMax<inf)
                subShape = SpectralLine(shapeParams);
            else
                subShape = HalfInfSpectralLine(shapeParams);
            end
            
            
            Int    = subShape.ComputeIntegrationVector();
            IP          = SubShapePts(this,subShape.Pts);
            weights = weight(subShape.Pts.y);
            weights = weights.';
            weights(isnan(weights) | weights == inf | weights == -inf) = 0;
            
            w = (Int .* weights) * IP;
            
            
        end
        
        
        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams)
            
            if(~isfield(shapeParams,'yMin'))
                shapeParams.yMin = 0;
            end
            if(~isfield(shapeParams,'yMax'))
                shapeParams.yMax = inf;
            end
            
            yMin = shapeParams.yMin;
            yMax = shapeParams.yMax;
            
            if(yMin>-inf && yMax<inf)
                geom.N = shapeParams.N;
                geom.yMin = yMin;
                geom.yMax = yMax;
                subShape = SpectralLine(geom);
            elseif(yMin>-inf && yMax == inf)
                geom.N = shapeParams.N;
                geom.yMin = abs(shapeParams.yMin);
                geom.L = abs(shapeParams.L);
                subShape = HalfInfSpectralLine(geom);
            elseif(yMin==-inf && yMax < inf)
                geom.N = shapeParams.N;
                geom.yMin = -abs(shapeParams.yMin);
                geom.L    = -abs(shapeParams.L);
                subShape = HalfInfSpectralLine(geom);
            else
                geom.N   = shapeParams.N;
                geom.L   = shapeParams.L;
                subShape = InfSpectralLine(geom);
            end
                
            
            % initialize matrix
            fP = f(0,subShape.Pts.y);

            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end

            nf = size(fP,2);
            M_conv = zeros(this.N,this.N,nf);  
                       
            % Integration and interpolation on subshape
            Int    = subShape.ComputeIntegrationVector();
            IntRep = Int(ones(1,nf),:);

            % interpolation onto subshape from the full domain
            IP          = SubShapePts(this,subShape.Pts);

            thisY = this.Pts.y;
            subY = subShape.Pts.y;
            
            for i=1:this.N

                % note that in the spherical case the functions depend on
                % both radial coordinates, rather than just the separation
                fP = f(thisY(i),subY);
                
                if(nDim>2)
                    fP_old = fP;
                    fP = reshape(fP,size(fP,1),[]);
                end
                
                IntFTemp = (IntRep.*fP.')*IP;

                if(size(IntFTemp,1)>1)
                    IntFTemp = IntFTemp.';
                end

                M_conv(i,:,:) = IntFTemp;

            end

            M_conv(isnan(M_conv)) = 0;

            if(nDim>2)
                nfTemp = size(fP_old);
                nfTemp = nfTemp(2:end);
                M_conv = reshape(M_conv,[ this.N, this.N, nfTemp ]);
            end
 
        end                        
        function M_conv = ComputeConvolutionMatrixPointwise(this,f,shapeParams)
           
            yMin = shapeParams.yMin;
            yMax = shapeParams.yMax;
            
            if(yMin>-inf && yMax<inf)
                geom.N = shapeParams.N;
                geom.yMin = yMin;
                geom.yMax = yMax;
                subShape = SpectralLine(geom);
            elseif(yMin>-inf && yMax == inf)
                geom.N = shapeParams.N;
                geom.yMin = abs(shapeParams.yMin);
                geom.L = abs(shapeParams.L);
                subShape = HalfInfSpectralLine(geom);
            elseif(yMin==-inf && yMax < inf)
                geom.N = shapeParams.N;
                geom.yMin = -abs(shapeParams.yMin);
                geom.L    = -abs(shapeParams.L);
                subShape = HalfInfSpectralLine(geom);
            else
                geom.N   = shapeParams.N;
                geom.L   = shapeParams.L;
                subShape = InfSpectralLine(geom);
            end
                
            
            % initialize matrix
            fP = f(0,subShape.Pts.y);

            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end

            nf = size(fP,2);
            M_conv = zeros(this.N,this.N,nf);  
                       
            % Integration on subshape
            Int    = subShape.ComputeIntegrationVector();
            IntRep = Int(ones(1,nf),:);

            thisY = this.Pts.y;
            subY = subShape.Pts.y;
                       
            for i=1:this.N

                shiftedSubShapeY = subY + thisY(i); 
                subShape.Pts.y = shiftedSubShapeY;
                % interpolation onto subshape from the full domain
                IP  = SubShapePts(this,subShape.Pts);
                
                % note that in the spherical case the functions depend on
                % both radial coordinates, rather than just the separation
                fP = f(thisY(i),shiftedSubShapeY);
                
                if(nDim>2)
                    fP_old = fP;
                    fP = reshape(fP,size(fP,1),[]);
                end
                
                IntFTemp = (IntRep.*fP.')*IP;

                if(size(IntFTemp,1)>1)
                    IntFTemp = IntFTemp.';
                end

                M_conv(i,:,:) = IntFTemp;

            end

            M_conv(isnan(M_conv)) = 0;

            if(nDim>2)
                nfTemp = size(fP_old);
                nfTemp = nfTemp(2:end);
                M_conv = reshape(M_conv,[ this.N, this.N, nfTemp ]);
            end
 
        end        
        
        function plot(this,V,options)
        
            if(options.dist)
                weight = diag(4*pi*this.Pts.y.^2);
            
                V = weight*V;
                V(isnan(V)) = 0;
            end
            
            plot@Interval(this,V,options)
            
        end
        
   end
   
end
