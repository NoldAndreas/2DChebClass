classdef (Abstract) Polar_SpectralFourierNoOrigin < SpectralFourier
    %y1 - radial variable r
    %y2 - angular variable theta
 
    methods 
        function this = Polar_SpectralFourierNoOrigin(N1,N2)
            this@SpectralFourier(N1,N2);                                    
                        
            this.polar       = 'polar';
        end
    end
    
    methods (Access = public)            
        function Int = ComputeIntegrationVector(this)
            
            Int = ComputeIntegrationVector@SpectralFourier(this);
            
            r = this.Pts.y1_kv';
            r(r==inf) = 0;
            Int(Int==inf) = 0;
            
            Int = Int.*r;

            this.Int = Int;
        end  
        function Diff = ComputeDifferentiationMatrix(this)            
            
            Diff = ComputeDifferentiationMatrix@SpectralFourier...
                                                (this,this.Pts);                           

            Diff   = PolarDiffOperators(Diff,this.Pts); 
            
            this.Diff = Diff;
        end                
        function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)            
            Interp = ComputeInterpolationMatrix@SpectralFourier(this,interp1,interp2,fullInterpGrid,saveBool);            
        end        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            N1 = this.N1;   N2 = this.N2;

            rs  = this.Pts.y1_kv; ts  = this.Pts.y2_kv;
            xs  = rs.*cos(ts);
            ys  = rs.*sin(ts);            

            M_conv = zeros(N1*N2,N1*N2);

            for i=1:(N1*N2)

                xd  = xs(i)-xs;
                yd  = ys(i)-ys;

                [t_d,r_d] = cart2pol(xd,yd);

                r_d((rs(i) == inf) | (rs ==inf)) = inf;
                r_d((rs(i) == inf) & (rs ==inf) & (t_d == 0)) = 0;            

                M_conv(i,:)  = this.Int.*f(r_d,t_d)';
            end   
            if((nargin== 3) && saveBool)
                this.Conv = M_conv;
            end
                
        end
        function Ind = ComputeIndices(this)

            x_kv    = kron(this.Pts.x1,ones(size(this.Pts.x2)));
            outR    = (x_kv == 1);        
            inR     = (x_kv == -1);        
            
            infBoundary = any(~isfinite(this.Pts.y1_kv));
            finite   = (~infBoundary & outR) | inR;
            finite1  = finite;
            finite2  = false(size(x_kv));
            infinite =  infBoundary & outR;
            
            Z                   = zeros(this.N1*this.N2);        
            nRout               = Z;
            nRin                = Z;

            nRout(outR,outR)    = speye(sum(outR));
            nRin(inR,inR)       = -speye(sum(inR));

            normalOutR = sparse( [nRout(outR,:)  Z(outR,:)]);
            normalInR  = sparse( [nRin(inR,:)  Z(inR,:)]);
            
            normalFinite   = sparse( [ (~infBoundary*nRout(finite,:) + nRin(finite,:) )  Z(finite,:)]);
            normalFinite1  = sparse( (~infBoundary*nRout(finite,:) + nRin(finite,:) ) );
            normalFinite2  = sparse( zeros(0,length(x_kv)) );
            
            normalInfinite = sparse( [ infBoundary*nRout(infinite,:)  Z(infinite,:)]);
            
            Ind = struct('outR',outR,...
                       'normalOutR',normalOutR ,...
                       'inR',inR,...
                       'normalInR', normalInR,...
                       'bound', (outR | inR),...
                       'corners',false(this.N1*this.N2,1), ...
                       'finite',finite,'infinite',infinite, ...
                       'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
                       'finite1',finite1,'finite2',finite2,...
                       'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);
                   
            this.Ind = Ind;
        end
    end
       
end