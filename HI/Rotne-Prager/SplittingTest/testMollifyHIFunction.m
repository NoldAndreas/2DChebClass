function [C_Disc,C_Annulus,C_Full]=testMollifyHIFunction

    clear all
    close all

    Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[10,10],'L1',4,...
                           'y2Min',-inf,'y2Max',inf,'L2',4);

    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                        'y2Min',-5,'y2Max',5,'N2',100);

    IDC = InfSpace_FMT(Phys_Area);

    [Pts,Diff,Int,Ind,~]  = IDC.ComputeAll(Plot_Area);

   
    RMax = 20;
    
    params.width = 1;
    
    params.N = [10;10];
    params.sigmaHS = 0.5;
    
    
    params.R = params.sigmaHS;
    
    C_Disc = IDC.ComputeConvolutionMatrix(@doMollify,params); 
    
     params = rmfield(params,'R');
     params.L = 2;
%    params.R = RMax;
    
    C_Full = IDC.ComputeConvolutionMatrix(@doMollify,params);
    
%     params = rmfield(params,'R');
%     params.RMax = RMax;
    params.RMin = params.sigmaHS;
    
    C_Annulus = IDC.ComputeConvolutionMatrix(@doMollify,params); 
    
    C_Test = C_Full - C_Disc;
    max(max(C_Test - C_Annulus))
    
    g = g1(IDC.Pts.y1_kv,IDC.Pts.y2_kv);
    
    C_Annulus = C_Annulus(:,:,1,1);
    C_Test = C_Test(:,:,1,1);
    
    %[C_Annulus*g, C_Test*g]
    
    [(C_Annulus*g - C_Test*g), sqrt(IDC.Pts.y1_kv.^2 + IDC.Pts.y2_kv.^2)]
    
    %----------------------------------------------------------------------
    
    function z = doMollify(x,y) 
        params.HIfn = str2func('HITestG');
        %params.HIfn = str2func('RP12_2D');
        z = MollifyHIFunction(x,y,params);
    end

    function g = g1(X,Y)
       sigmaX = 1; sigmaY = 1;
       X0 = 1;  Y0 = 1;
       g = ( 2*pi *sigmaX*sigmaY ).^(-1) * exp( - ( (X-X0).^2/2/sigmaX^2 + (Y-Y0).^2/2/sigmaY^2 ) );
       g(isnan(X) | isnan(Y)) = 0;
    end

end