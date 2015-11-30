function ThesisNumerics_NumericsPaper_GaussianInterpolation    

    AddPaths('CodePaper');
    set(0,'defaultaxesfontsize',20);
    set(0,'defaultlinelinewidth',2);

    defaultPos = [0 0 1000 800];

    saveDir = [pwd filesep 'Computations' filesep 'CodePaper' filesep 'Images' filesep];

    close all;

    geom.N = 50;
    geom.L = 2;  % optimize by looking at a derivative, or do parameter testing
    
    infLine = InfSpectralLine(geom);
    
    % for Gaussian plots
    PlotArea = struct('N',1500,'yMin',-25,'yMax',25);
    
    % for kernel plots
    %PlotArea = struct('N',1000,'yMin',-10,'yMax',10);
    
    [Pts,Diff,Int,Ind,Interp] = infLine.ComputeAll(PlotArea);

    % this case breaks the basic convolution
    mu1 = 0;
    mu2 = 2;   
    sigma1 = 0.25;
    sigma2 = 4;

    mu3 = 1;
    sigma3 = 1;
    
    % this case is pretty good for both (still better for ptwise
%     mu1 = 0;
%     mu2 = 1;    
%     sigma1 = 0.25;
%     sigma2 = 1;


    muConv12 = mu1+mu2;
    sigmaConv12 = sqrt(sigma1^2 + sigma2^2);

    muConv13 = mu1+mu3;
    sigmaConv13 = sqrt(sigma1^2 + sigma3^2);

    y = Pts.y;
    sigma1 = 1;
    
    figure('Position',[0 0 450 200],'color','white');  
    
    subplot(1,3,1);
    k = geom.N-20;
    yTest = y(k);
    shiftedKernel = Gaussian1(y-yTest);
    plot(Interp.pts,Interp.InterPol*shiftedKernel,'k'); hold on;
    plot(y,shiftedKernel,'ok');
    xlim([-25,25]);  ylim([-0.1,0.5]);
    xlabel('$y$','interpreter','latex');
    title(['$g(y-y_{',num2str(k),'})$'],'Interpreter','latex');
    pbaspect([1 1 1]);
    
    subplot(1,3,2);
    k = geom.N-6;
    yTest = y(k);    
    shiftedKernel = Gaussian1(y-yTest);
    plot(Interp.pts,Interp.InterPol*shiftedKernel,'k');hold on;
    plot(y,shiftedKernel,'ok');
    xlim([-25,25]);  ylim([-0.1,0.5]);
    xlabel('$y$','interpreter','latex');
    title(['$g(y-y_{',num2str(k),'})$'],'interpreter','latex');
    pbaspect([1 1 1]);
    
    subplot(1,3,3);
    k = geom.N-3;
    yTest = y(k);   
    shiftedKernel = Gaussian1(y-yTest);
    plot(Interp.pts,Interp.InterPol*shiftedKernel,'k'); hold on;
    plot(y,shiftedKernel,'ok');
    pbaspect([1 1 1]);
    
    xlim([-25,25]);  ylim([-0.1,0.5]);
    xlabel('$y$','interpreter','latex');
    title(['$g(y-y_{',num2str(k),'})$'],'interpreter','latex');
    SaveFigure('kernel_all');
    
    sigma1 = 0.25;
    g1 = Gaussian1(y);
    g2 = Gaussian2(y);
    g3 = Gaussian3(y);
    
    gConv12 = GaussianConv12(y);
    gConv13 = GaussianConv13(y);
    
    shapeParams.N  = 50;
    shapeParams.L  = 2;
    
    % false gives pointwise convolution
    convMatrixPtwise = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,false);
    
    % true uses naive convolution
    convMatrixBasic = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,true);

    gPtwise2 = convMatrixPtwise*g2;
    gBasic2 = convMatrixBasic*g2;

    gPtwise3 = convMatrixPtwise*g3;
    gBasic3 = convMatrixBasic*g3;

     
%     figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*g1,'r');
%     hold on
%     plot(Interp.pts,Interp.InterPol*g2,'b');
%     plot(Interp.pts,Interp.InterPol*g3,'m');
%     plot(y,g1,'or');
%     plot(y,g2,'ob');
%     plot(y,g3,'om');
%     xlim([-20,20]);
%     xlabel('$y$','interpreter','latex');
%     ylabel('$g_j$','interpreter','latex');
%     save2pdf([saveDir 'convolutionGaussians.pdf'],gcf);
%     
%     figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*gPtwise2,'b');
%     hold on
%     plot(y,gPtwise2,'ob');
%     plot(Interp.pts,Interp.InterPol*gBasic2,'--b');
%     plot(y,gBasic2,'xb');
% 
%     plot(Interp.pts,Interp.InterPol*gPtwise3,'m');
%     hold on
%     plot(y,gPtwise3,'om');
%     plot(Interp.pts,Interp.InterPol*gBasic3,'--m');
%     plot(y,gBasic3,'xm');
% 
%     
%     %plot(y,gConv,'og');
%     xlim([-20,20]);
%     xlabel('$y$','interpreter','latex');
%     ylabel('$g_1 \ast g_2$','interpreter','latex');
%     save2pdf([saveDir 'convolutionResults.pdf'],gcf);
    
    ptwiseErr2 = L2norm(gPtwise2-gConv12)/L2norm(gConv12)
    basicErr2  = L2norm(gBasic2-gConv12)/L2norm(gConv12)
        
    ptwiseErr3 = L2norm(gPtwise3-gConv13)/L2norm(gConv13)
    basicErr3  = L2norm(gBasic3-gConv13)/L2norm(gConv13)

    
    function g = Gaussian1(x)
        g = 1/sqrt(2*pi)/sigma1 * exp( -(x-mu1).^2/(2*sigma1^2) );
    end

    function g = Gaussian2(x)
        g = 1/sqrt(2*pi)/sigma2 * exp( -(x-mu2).^2/(2*sigma2^2) );
    end

    function g = Gaussian3(x)
        g = 1/sqrt(2*pi)/sigma3 * exp( -(x-mu3).^2/(2*sigma3^2) );
    end


    function g = GaussianConv12(x)
        g = 1/sqrt(2*pi)/sigmaConv12 * exp( -(x-muConv12).^2/(2*sigmaConv12^2) );
    end

    function g = GaussianConv13(x)
        g = 1/sqrt(2*pi)/sigmaConv13 * exp( -(x-muConv13).^2/(2*sigmaConv13^2) );
    end


    function Norm = L2norm(f)
        Norm = Int*(f.^2);
    end


end