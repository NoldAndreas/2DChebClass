        function testAll(this)
            
            int = this.ComputeIntegrationVector;
            
            diff = this.ComputeDifferentiationMatrix;
            
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;
            
            g = 1/pi * exp( - (y1.^2 + y2.^2) );
            
            intError = (int*g -1)
            
            DError = max(abs(diff.Dy1*g - (-g.*2.*y1)))
           
            interp1 = (-1:0.01:1).';
            interp2 = interp1;

            Interp = this.ComputeInterpolationMatrix(interp1,interp2,true);
            
            gInterp = Interp.InterPol*g;
            
            gPlot =  1/pi * exp( - (Interp.pts1.^2 + Interp.pts2.^2) );
            
            interpError = max( abs ( gInterp - gPlot) )
            
            
            sigma1 = 1;  sigma2  = 2;
            mu1    = 0;  mu2     = 5;
            
            r = sqrt( y1.^2 + y2.^2 );
            
            M_conv = ComputeConvolutionMatrix(this,@g1);
            
            conv = M_conv*g2(y1,y2);
            
            test = exp( - ( r - (mu1+mu2) ).^2 / ( 2 * (sigma1^2 + sigma2^2) ) );
            test = test * 1 / ( 2*pi * (sigma1^2 + sigma2^2) );
            
            convRelError = max(abs(conv-test)./abs(test))
            
            temp = reshape(conv - test,this.N1,this.N2);
            surf(this.Pts.x1,this.Pts.x2,abs(temp));
             
%             surf(this.Pts.x1,this.Pts.x2,reshape(conv,this.N1,this.N2),'FaceColor','r');
%             hold on
%             surf(this.Pts.x1,this.Pts.x2,reshape(test,this.N1,this.N2),'FaceColor','b');

            
            % aux functions
            function g = g1(x,y)
               sigma = sigma1;
               mu    = mu1;
               
               R = sqrt(x.^2 + y.^2);
               
               g = ( 2*pi *sigma^2 ).^(-1) * exp( - (R-mu).^2/2/sigma^2 );
            end

            function g = g2(x,y)
               sigma = sigma2;
               mu    = mu2;
               
               R = sqrt(x.^2 + y.^2);
               
               g = ( 2*pi *sigma^2 ).^(-1) * exp( - (R-mu).^2/2/sigma^2 );
            end
            
        end
        