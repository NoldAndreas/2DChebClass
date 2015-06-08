classdef (Abstract) Polar_SpectralFourier < SpectralFourier
    %y1 - radial variable r
    %y2 - angular variable theta
    
    methods 
        function this = Polar_SpectralFourier(N1,N2)
            this@SpectralFourier(2*N1,N2);                                    
            
            this.N1          = N1;
            this.Pts.x1full  = this.Pts.x1;
            this.Pts.x1      = this.Pts.x1(N1+1:end);  
            
            this.Pts.N1 = N1; this.Pts.N2 = N2;
                        
            this.FFTMatrix   = this.FFTMatrix(1:end/2,1:end/2);
            this.IFFTMatrix  = this.IFFTMatrix(1:end/2,1:end/2);   
            
            this.polar       = 'polar';                        
        end
    end
    
    methods (Access = public)            
       function Int = ComputeIntegrationVector(this,t1Odd)
            
            N1 = this.N1;
            N2 = this.N2;
            %for values at infinity: zero required!
            %assumes that convergence is good enough for integration and
            %interpolation

            [h1,dy2] = this.PhysSpace2(this.Pts.x2);
            %1st step: integrate in phi-direction. This gives us values for
            %positive r
            %M1 = (dy2').*conj(fft(ones(1,N2)))/(N2^2);
            M1 = (dy2').*ones(1,N2)/N2;

            %2nd step: obtain values for full r-range, and multiply with r
            y1 = this.PhysSpace1(this.Pts.x1);
            M2 = [fliplr(diag(-y1))' ; diag(y1)];
            %x1full = [-flipud(this.Pts.x1) ; this.Pts.x1];

            %3rd step: Interpolate onto a Chebychev grid only in the positive
            %domain, having the same no of points,
            [xPos,wIntPos] = ClenCurtFlip(N1-1);     %Version 1
            
            if((nargin == 2) && t1Odd)
                wInt1O  = pi*[0.5 ones(1,N1-2) 0.5]/(N1-1);%sin(theta)!! 
                t       = pi*(0:(N1-1))/(N1-1);
                wIntPos   = wInt1O.*sin(t);
            end
            %[xPos,wIntPos] = ClenCurtFlip(4*(N1-1));   %Version 2

            %[y1New,dy1New] = Maps.PhysSpace1((1+xPos)/2);  %Version 1
            %dy1New         = dy1New/2;  %Version 1
            %x1New          = (1+xPos)/2;  %Version 1

            L05_1New          = this.PhysSpace1(sqrt(0.5)); 
            if(this.PhysSpace1(1) == inf)
                [y1New,dy1New]    = QuotientMap(xPos,L05_1New,0,inf);
            else
                [y1New,dy1New]    = LinearMap(xPos,0,this.PhysSpace1(1));%QuadraticMapRight(xPos);%%%LinearMap(xPos,0,this.PhysSpace1(1));% % QuadraticMap(xPos,this.PhysSpace1(1)); %QuotientMap(xPos,L05_1New,0,this.PhysSpace1(1));%LinearMap(xPos,0,this.PhysSpace1(1));%; %Version 2 % 
            end
            %y1New = y1New*this.PhysSpace1(1); dy1New = dy1New*this.PhysSpace1(1);
            x1New             = this.CompSpace1(y1New);%Version 2
            M3                = barychebevalMatrix(this.Pts.x1full,x1New); 

            %4th step: Integrate on this new Chebychev Grid

            %5th step: take care of values at infinity
            dy1New(dy1New == inf)  = 0;
            M2(M2 == inf)          = 0;
            M2(M2 == -inf)         = 0;

            %6th step: combine steps in r- and in phi-direction        
            Int      = kron( (dy1New'.*wIntPos) *  M3 * M2 ,M1);
            this.Int = Int;
        end
       function Ind = ComputeIndices(this)

            x1_kv    = kron(this.Pts.x1,ones(size(this.Pts.x2)));
            outR    = (x1_kv == 1);        

            infBoundary = any(~isfinite(this.Pts.y1_kv));
            finite   = ~infBoundary & outR;
            finite1  = finite;
            finite2  = false(size(x1_kv));
            infinite =  infBoundary & outR;
            
            Z                   = zeros(this.N1*this.N2);        
            nRout               = Z;

            nRout(outR,outR)    = speye(sum(outR));
            
            normalOutR  = sparse( [nRout(outR,:)  Z(outR,:)]);
            
            normalFinite   = sparse( [~infBoundary*nRout(finite,:)  Z(finite,:)]);
            normalFinite1  = sparse( ~infBoundary*nRout(finite,:) );
            normalFinite2  = sparse( zeros(0,length(x1_kv)) );
            
            normalInfinite = sparse( [ infBoundary*nRout(infinite,:)  Z(infinite,:)]);
            
            Ind = struct('outR',outR,...
                       'normalOutR',normalOutR ,...
                       'bound',outR,...
                       'corners',false(this.N1*this.N2,1), ...
                       'finite',finite,'infinite',infinite, ...
                       'normalFinite',normalFinite,'normalInfinite',normalInfinite,...
                       'finite1',finite1,'finite2',finite2,...
                       'normalFinite1',normalFinite1,'normalFinite2',normalFinite2);
                   
            this.Ind = Ind;
        end
       function Diff = ComputeDifferentiationMatrix(this)            
            pts    = this.Pts;
            pts.x1 = this.Pts.x1full;
            mMask  = 1:(this.N1*this.N2);
            pMask  = (this.N1*this.N2+1):(2*this.N1*this.N2);
            
            Diff = ComputeDifferentiationMatrix@SpectralFourier...
                                                (this,pts);
            
            M           = kron(fliplr(eye(this.N1)),...
                                    circshift(eye(this.N2),[0,this.N2/2]));        

            Dy2Polar    = Diff.Dy2(pMask,pMask);
            DDy2Polar   = Diff.DDy2(pMask,pMask);
            Dy1Polar    = Diff.Dy1(pMask,mMask)*M + Diff.Dy1(pMask,pMask);        
            DDy1Polar   = Diff.DDy1(pMask,mMask)*M + Diff.DDy1(pMask,pMask);        
            Dy1Dy2Polar = Diff.Dy1Dy2(pMask,mMask)*M + Diff.Dy1Dy2(pMask,pMask);        

            DiffPolar_h = struct('Dy1',Dy1Polar,'Dy2',Dy2Polar,...
                                    'DDy1',DDy1Polar,'DDy2',DDy2Polar,...
                                    'Dy1Dy2',Dy1Dy2Polar);                                        

            Diff   = PolarDiffOperators(DiffPolar_h,this.Pts);            
            this.Diff = Diff;
        end                        
       function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)

            N1 = this.N1;       N2 = this.N2;
            
            Interp1     = barychebevalMatrix(this.Pts.x1full,interp1);

            %2nd variable: Fourier
            LL    = 2*pi;
            M     = floor((N2+1)/2);

            %%FFT Interpolation
            n21   = 1:(M-1);       % first half of points, excluding zero
            n22   = (M+1):(N2-1);  % second half

            % This uses that NT+1 is odd (ie NT is even)
            %calculate the values for positive r...
            phi       = interp2*2*pi;
            Interp2P  = [exp(2*pi*1i*phi/LL*0) , exp(2*pi*1i*phi/LL*n21) , cos(2*pi*M*phi/LL) , exp(2*pi*1i*phi/LL*(n22-N2))]/N2;

            %... and for negative r.
            phi      = mod(interp2*2*pi+pi,2*pi);
            Interp2N  = [exp(2*pi*1i*phi/LL*0) , exp(2*pi*1i*phi/LL*n21) , cos(2*pi*M*phi/LL) , exp(2*pi*1i*phi/LL*(n22-N2))]/N2;

            % kron form of the two interpolations                                                            

            Interpol1P = Interp1(:,N1+1:end);
            Interpol1N = Interp1(:,N1:-1:1);

            InterPolPolarOR = kron(Interpol1N,Interp2N) + kron(Interpol1P,Interp2P); 
            InterPolPolar = real(InterPolPolarOR*this.FFTMatrix);

            interp_y1 = this.PhysSpace1(interp1);
            interp_y2 = this.PhysSpace2(interp2);
            % kron form of the two interpolations                                                            
            Interp = struct('InterPol',InterPolPolar,...
                            'interp1',interp_y1,'interp2',interp_y2,...             
                            'Nplot1',length(interp1),...
                            'Nplot2',length(interp2),...
                            'N1',N1,'N2',N2);        
                        
            if((nargin >= 4) && fullInterpGrid)
                Interp.pts1 = kron(interp_y1,ones(size(interp_y2)));
                Interp.pts2 = kron(ones(size(interp_y1)),interp_y2);             
            end
            if((nargin >= 5) && saveBool)
                this.Interp = Interp;
            end
        end            
       function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            N1  = this.N1;  N2  = this.N2;
            Int = this.Int;  % 1 x N1*N2

            rs  = this.Pts.y1_kv; ts  = this.Pts.y2_kv;
            xs  = rs.*cos(ts);
            ys  = rs.*sin(ts);  
            
            fPTemp = f(rs);%,ts);
            
            fDim = size(fPTemp);
            
            nElts = prod(fDim(2:end));
            
            IntT = Int.';  % N1*N2 x 1
            
            IntT = IntT(:,ones(1,nElts)); % N1*N2 x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([N1*N2,N1*N2,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:(N1*N2) 
                
                xd  = xs(i)-xs;
                yd  = ys(i)-ys;

                [t_d,r_d] = cart2pol(xd,yd);

                r_d((rs(i) == inf) | (rs ==inf)) = inf;
                r_d((rs(i) == inf) & (rs ==inf) & (t_d == 0)) = 0;   
                
                fP          = f(r_d);%,t_d);
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end             
    end
       
end