 function [M_conv,M_conv_New] = ComputeConvolutionMatrix_Test(this,shapeParams)
            
            r = this.Pts.y1_kv;
            t = this.Pts.y2_kv;

            sigmax1 = 0.5 ;  sigmay1  = 1;
            sigmax2 = 10;  sigmay2  = 15;
            sigmax3 = 1;  sigmay3  = 0.7;

            x01 = 0 ; y01 = -2;
            x02 = 1  ; y02 = 0; 
            x03 = 0  ; y03 = 0.5; 
            
            function g = g1(r,t)
               sigmaX = sigmax1; sigmaY = sigmay1;
               X0 = x01;  Y0 = y01;
               X = r.*cos(t); Y = r.*sin(t);
               g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
               g(isnan(X) | isnan(Y)) = 0;
            end

            function g = g2(r,t)
               sigmaX = sigmax2; sigmaY = sigmay2;
               X0 = x02;  Y0 = y02;
               X = r.*cos(t); Y = r.*sin(t);
               g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
               g(isnan(X) | isnan(Y)) = 0;
            end
            
            function g = g3(r,t)
               sigmaX = sigmax3; sigmaY = sigmay3;
               X0 = x03;  Y0 = y03;
               X = r.*cos(t); Y = r.*sin(t);
               g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
               g(isnan(X) | isnan(Y)) = 0;
            end
            
            
            % infinite disc, convolution of two Gaussians
            shapeParams.N=[20;20]; shapeParams.L=2; 
            this.Int = this.ComputeIntegrationVector;
            
            M_conv     = ComputeConvolutionMatrix(this,@g1,[],true);
            M_conv_New = ComputeConvolutionMatrix(this,@g1,shapeParams);
            
            conv2 = M_conv*g2(r,t);
            conv_New2 = M_conv_New*g2(r,t);
            
            conv3 = M_conv*g3(r,t);
            conv_New3 = M_conv_New*g3(r,t);
            
            x = r.*cos(t); y = r.*sin(t);
 
            muX2 = x01+x02;  sigmaX2 = sqrt(sigmax1^2 + sigmax2^2);  
            muY2 = y01+y02;  sigmaY2 = sqrt(sigmay1^2 + sigmay2^2);
            exact2 = exp( - ( x - muX2 ).^2 / ( 2 * sigmaX2^2 ) ) ...
                    .* exp( - ( y - muY2 ).^2 / ( 2 * sigmaY2^2 ) );
            exact2 = exact2 * 1 / ( 2*pi * sigmaX2 * sigmaY2 );
            exact2(isnan(exact2))=0;
            
            muX3 = x01+x03;  sigmaX3 = sqrt(sigmax1^2 + sigmax3^2);  
            muY3 = y01+y03;  sigmaY3 = sqrt(sigmay1^2 + sigmay3^2);
            exact3 = exp( - ( x - muX3 ).^2 / ( 2 * sigmaX3^2 ) ) ...
                    .* exp( - ( y - muY3 ).^2 / ( 2 * sigmaY3^2 ) );
            exact3 = exact3 * 1 / ( 2*pi * sigmaX3 * sigmaY3 );
            exact3(isnan(exact3))=0;
            
            rRange = (-0.9:0.05:0.9)';
            Interp   = this.ComputeInterpolationMatrix(rRange,(0:0.02:1)',true,false);
            
            figure
            doPlots(this,exact2);
            set(gcf,'Name','Exact 2');
            
            figure
            doPlots(this,conv2);
            set(gcf,'Name','Standard Convolution 2');
            
            figure
            doPlots(this,conv_New2);
            set(gcf,'Name','Pointwise Convolution 2');
            
            disp(['InfDisc: ComputeConvolutionMatrix_Test ''standard'' error: ' ...
                      num2str( sum(abs(exact2-conv2).^2) / sum(abs(exact2).^2) )])
            disp(['InfDisc: ComputeConvolutionMatrix_Test ''pointwise'' error: ' ...
                      num2str( sum(abs(exact2-conv_New2).^2) / sum(abs(exact2).^2) )])
            
            figure
            doPlots(this,exact3);
            set(gcf,'Name','Exact 3');
            
            figure
            doPlots(this,conv3);
            set(gcf,'Name','Standard Convolution 3');
            
            figure
            doPlots(this,conv_New3);
            set(gcf,'Name','Pointwise Convolution 3');
            
            disp(['InfDisc: ComputeConvolutionMatrix_Test ''standard'' error: ' ...
                      num2str( sum(abs(exact3-conv3).^2) / sum(abs(exact3).^2) )])
            disp(['InfDisc: ComputeConvolutionMatrix_Test ''pointwise'' error: ' ...
                      num2str( sum(abs(exact3-conv_New3).^2) / sum(abs(exact3).^2) )])
            
        end