classdef Sphere_Old < Polar_SpectralFourier

    properties        
        R        
        Ch         
    end
    
    methods        
        function this = Sphere_old(Geometry)
            this@Polar_SpectralFourier(Geometry.N(1),Geometry.N(2));
            
            this.R      = Geometry.R;             
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            InitializationPts(this);                                    
        end                         
    end
    
    methods (Access = public)
        function [r,dr] = PhysSpace1(this,x)%,dx,ddx,dddx,ddddx]
            %[r,dr,dx,ddx,dddx,ddddx] = LinearMap(x,-this.R,this.R);
            [r,dr] = QuadraticMap(x);
            r      = this.R*r;
            dr     = this.R*dr;
        end
        function xf = CompSpace1(this,r)
            r   = r/this.R;            
            xf = bisection(@QuadraticMap,-1,1,r);%,options)            
            %xf  = fzero(fun,r);
            %xf = InvLinearMap(r,-this.R,this.R);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,xT)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end
        function xf = CompSpace2(this,th)
            xf = th/(2*pi);
        end                        
        
        function [Int,area] = ComputeIntegrationVector(this,t1Odd)
            
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
            
            if(((nargin == 2) && t1Odd))
                wInt1O  = pi*[0.5 ones(1,N1-2) 0.5]/(N1-1);%sin(theta)!! 
                t       = pi*(0:(N1-1))/(N1-1);
                wIntPos   = wInt1O.*sin(t);
            end
            %[xPos,wIntPos] = ClenCurtFlip(4*(N1-1));   %Version 2

            %[y1New,dy1New] = Maps.PhysSpace1((1+xPos)/2);  %Version 1
            %dy1New         = dy1New/2;  %Version 1
            %x1New          = (1+xPos)/2;  %Version 1

            L05_1New          = this.PhysSpace1(sqrt(0.5)); 
            [y1New,dy1New]    = QuadraticMapRight(xPos);%LinearMap(xPos,0,this.PhysSpace1(1));%QuadraticMapRight(xPos);%LinearMap(xPos,0,this.PhysSpace1(1));%QuadraticMapRight(xPos);%%%LinearMap(xPos,0,this.PhysSpace1(1));% % QuadraticMap(xPos,this.PhysSpace1(1)); %QuotientMap(xPos,L05_1New,0,this.PhysSpace1(1));%LinearMap(xPos,0,this.PhysSpace1(1));%; %Version 2 % 
            y1New = y1New*this.PhysSpace1(1); dy1New = dy1New*this.PhysSpace1(1);
            x1New             = this.CompSpace1(y1New);%Version 2
            M3                = barychebevalMatrix(this.Pts.x1full,x1New); 

            %4th step: Integrate on this new Chebychev Grid

            %5th step: take care of values at infinity
            dy1New(dy1New == inf)  = 0;
            M2(M2 == inf)          = 0;
            M2(M2 == -inf)         = 0;

            %6th step: combine steps in r- and in phi-direction        
            Int      = kron( (dy1New'.*wIntPos) *  M3 * M2 ,M1);
            
            Int      = 2*Int.*sqrt(this.R^2 - this.Pts.y1_kv.^2)';
            if(sum(imag(Int)) > 0)
                error('Disc:ComputeIntegrationVector: Imaginary part in integration vector');
            end
            this.Int = Int;
            area     = 4/3*pi*this.R^3;
            
             if(nargout < 2)
                disp(['Sphere: Error of integration of area (ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);                
            end
            
        end
    end
end