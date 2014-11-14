function ComputeConvolutionMatrix_Test_Short_2(this)

    % Run this with 
    %    infGeom.N=[20;20]; infGeom.L1 = 10; infGeom.L2 = 10;
    % and   sigmax4 = 10;  sigmay4  = 10;
    % to see the difference between the standard and
    % pointwise convolution.
    %
    % Or run with 
    %    infGeom.N=[20;20]; infGeom.L1 = 2; infGeom.L2 = 2;
    % and   sigmax4 = 1;  sigmay4  = 1;
    % to see the accuracy of both approaches
   
    X = this.Pts.y1_kv;
    Y = this.Pts.y2_kv;

    sigmax1 = 1 ;  sigmay1  = 1;
    sigmax2 = 0.5;  sigmay2  = 0.5;
    sigmax3 = 1;  sigmay3  = 1;
    %sigmax4 = 1;  sigmay4  = 1;
    sigmax4 = 10;  sigmay4  = 10;

    x01 = -1 ; y01 = 0;
    x02 = 0  ; y02 = 0; 
    x03 = 0  ; y03 = 0; 
    x04 = 1  ; y04 = 0.5; 

    muX14 = x01+x04;  sigmaX14 = sqrt(sigmax1^2 + sigmax4^2);  
    muY14 = y01+y04;  sigmaY14 = sqrt(sigmay1^2 + sigmay4^2);
    exact14 = exp( - ( X - muX14 ).^2 / ( 2 * sigmaX14^2 ) ) ...
            .* exp( - ( Y - muY14 ).^2 / ( 2 * sigmaY14^2 ) );
    exact14 = exact14 * 1 / ( 2*pi * sigmaX14 * sigmaY14 );
    exact14(isnan(exact14))=0;

    muX24 = x02+x04;  sigmaX24 = sqrt(sigmax2^2 + sigmax4^2);  
    muY24 = y02+y04;  sigmaY24 = sqrt(sigmay2^2 + sigmay4^2);
    exact24 = exp( - ( X - muX24 ).^2 / ( 2 * sigmaX24^2 ) ) ...
            .* exp( - ( Y - muY24 ).^2 / ( 2 * sigmaY24^2 ) );
    exact24 = exact24 * 1 / ( 2*pi * sigmaX24 * sigmaY24 );
    exact24(isnan(exact24))=0;

    muX34 = x03+x04;  sigmaX34 = sqrt(sigmax3^2 + sigmax4^2);  
    muY34 = y03+y04;  sigmaY34 = sqrt(sigmay3^2 + sigmay4^2);
    exact34 = exp( - ( X - muX34 ).^2 / ( 2 * sigmaX34^2 ) ) ...
            .* exp( - ( Y - muY34 ).^2 / ( 2 * sigmaY34^2 ) );
    exact34 = exact34 * 1 / ( 2*pi * sigmaX34 * sigmaY34 );
    exact34(isnan(exact34))=0;

    XRange = (-0.9:0.05:0.9)';
    YRange = (-0.9:0.05:0.9)';
    this.ComputeInterpolationMatrix(XRange,YRange,true,true);
   
    % infinite disc, convolution of two Gaussians
    shapeParams.N=[20;20]; shapeParams.L=2; 
    this.Int = this.ComputeIntegrationVector;

    M_conv_M = ComputeConvolutionMatrix(this,@gMat,[],true);

    M_conv_New_M = ComputeConvolutionMatrix(this,@gMat,shapeParams);

    
    R = 0.5;
    discParams.N = [20;20];
    discParams.R = R;
    
    M_conv_inner_M =  ComputeConvolutionMatrix(this,@gMat,discParams);

    annulusParams.N = [20;20];
    annulusParams.L = 2;
    annulusParams.RMin = R;
    
    M_conv_outer_M =  ComputeConvolutionMatrix(this,@gMat,annulusParams);
    
    M_conv_New_M_2 = M_conv_inner_M + M_conv_outer_M;
    
    
    N1 = this.Pts.N1;
    N2 = this.Pts.N2;

    conv_M = zeros(N1*N2,2,2);
    conv_New_M = zeros(N1*N2,2,2);
    conv_New_M_2 = zeros(N1*N2,2,2);

    for k1 = 1:2
        for k2 = 1:2 
            conv_M(:,k1,k2)     = M_conv_M(:,:,k1,k2)*g4(X,Y);
            conv_New_M(:,k1,k2) = M_conv_New_M(:,:,k1,k2)*g4(X,Y);
            conv_New_M_2(:,k1,k2) = M_conv_New_M_2(:,:,k1,k2)*g4(X,Y);
        end
    end


    figure
    this.plot([exact14 exact24 exact34]);
    set(gcf,'Name','Exact');

    figure
    this.plot(reshape(conv_New_M,N1*N2,4));
    set(gcf,'Name','Pointwise Convolution Matrix');
    
    figure
    this.plot(reshape(conv_New_M_2,N1*N2,4));
    set(gcf,'Name','Split Convolution Matrix');

    figure
    this.plot(reshape(conv_M,N1*N2,4));
    set(gcf,'Name','Standard Convolution Matrix');

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''pointwise'' error: ' ...
              num2str( sum(abs(exact14-conv_New_M(:,1,1)).^2) / sum(abs(exact14).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''pointwise'' error: ' ...
            num2str( sum(abs(exact24-conv_New_M(:,1,2)).^2) / sum(abs(exact24).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''pointwise'' error: ' ...
            num2str( sum(abs(exact24-conv_New_M(:,2,1)).^2) / sum(abs(exact24).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''pointwise'' error: ' ...
            num2str( sum(abs(exact34-conv_New_M(:,2,2)).^2) / sum(abs(exact34).^2) )])

        
    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''split'' error: ' ...
              num2str( sum(abs(exact14-conv_New_M_2(:,1,1)).^2) / sum(abs(exact14).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''split'' error: ' ...
            num2str( sum(abs(exact24-conv_New_M_2(:,1,2)).^2) / sum(abs(exact24).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''split'' error: ' ...
            num2str( sum(abs(exact24-conv_New_M_2(:,2,1)).^2) / sum(abs(exact24).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''split'' error: ' ...
            num2str( sum(abs(exact34-conv_New_M_2(:,2,2)).^2) / sum(abs(exact34).^2) )])
        
        
    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''standard'' error: ' ...
              num2str( sum(abs(exact14-conv_M(:,1,1)).^2) / sum(abs(exact14).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''standard'' error: ' ...
            num2str( sum(abs(exact24-conv_M(:,1,2)).^2) / sum(abs(exact24).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''standard'' error: ' ...
            num2str( sum(abs(exact24-conv_M(:,2,1)).^2) / sum(abs(exact24).^2) )])

    disp(['InfSpace: ComputeConvolutionMatrix_Test Matrix ''standard'' error: ' ...
            num2str( sum(abs(exact34-conv_M(:,2,2)).^2) / sum(abs(exact34).^2) )])

%--------------------------------------------------------------------------        
        
    function g = g1g2(X,Y)
        g = [g1(X,Y) g2(X,Y)];
    end

    function g = gMat(X,Y)
        g = zeros(length(X),2,2);
        g(:,1,1) = g1(X,Y);
        g(:,1,2) = g2(X,Y);
        g(:,2,1) = g2(X,Y);
        g(:,2,2) = g3(X,Y);
    end


    function g = g1(X,Y)
       sigmaX = sigmax1; sigmaY = sigmay1;
       X0 = x01;  Y0 = y01;
       g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
       g(isnan(X) | isnan(Y)) = 0;
    end

    function g = g2(X,Y)
       sigmaX = sigmax2; sigmaY = sigmay2;
       X0 = x02;  Y0 = y02;
       g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
       g(isnan(X) | isnan(Y)) = 0;
    end

    function g = g3(X,Y)
       sigmaX = sigmax3; sigmaY = sigmay3;
       X0 = x03;  Y0 = y03;
       g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
       g(isnan(X) | isnan(Y)) = 0;
    end

    function g = g4(X,Y)
       sigmaX = sigmax4; sigmaY = sigmay4;
       X0 = x04;  Y0 = y04;
       g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
       g(isnan(X) | isnan(Y)) = 0;
    end
        
        
end
      