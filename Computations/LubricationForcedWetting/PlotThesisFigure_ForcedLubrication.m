function PlotThesisFigure_ForcedLubrication(y,hP,lambda,delta,thetaAP,Cout,Cin,h0P,h1P,h2P,hIP_0,hIP_1,hIP_2,HIS,yInner,ySInner,IP_In,IP_SIn)



        figure('color','white','Position',[0 0 400 400]);
      
        lw_l = 1;        
        yP   = 10.^((-5.6:0.01:1)');   
        
        yP_O = yP(yP>1e-4);
        plot(yP_O,thetaAP+(delta/thetaAP^2)*(log(yP_O)+Cout),'r-.','linewidth',1.); hold on;
        plot(yP_O,thetaAP+(delta/thetaAP^2)*(log(yP_O)+Cout)+...
                        (delta^2/thetaAP^5)*(-(log(yP_O)).^2-2*Cout*log(yP_O)),...
                        'r:','linewidth',1.);  hold on;
                  
        maskY =(y>1e-4);
        plot(y(maskY),h0P(maskY),'b--','linewidth',lw_l); hold on;
        plot(y(maskY),h0P(maskY)+delta*h1P(maskY),'b-.','linewidth',lw_l);    
        plot(y(maskY),h0P(maskY)+delta*h1P(maskY)+delta^2*h2P(maskY),'b:','linewidth',lw_l);                

        mark_h  = (y*lambda<1e-4);
        mark_h4 = (y*lambda<1e-4);
        plot(y(mark_h4)*lambda,ones(size(y(mark_h4))),'b--','linewidth',lw_l);                
        plot(y(mark_h4)*lambda,hIP_0(mark_h4)+delta*hIP_1(mark_h4),'b-.','linewidth',lw_l);
        plot(y(mark_h)*lambda,hIP_0(mark_h) + ...
                               delta*hIP_1(mark_h) + ...
                               delta^2*hIP_2(mark_h),'b:','linewidth',lw_l);        

        HIS.plot(hP,struct('plain',true,'linecolor','k','linewidth',lw_l));                           
        plot(yInner,IP_In.InterPol*hP,'k','linewidth',lw_l); hold on;
        plot(ySInner,IP_SIn.InterPol*hP,'k','linewidth',lw_l);

        %Intermediate Region result
        plot(yP,(thetaAP^3+3*delta*(log(yP)+Cout)).^(1/3),'g--','linewidth',1.);    
        
        %plot(lambda*res.L100.y,res.L100.hP,'k','linewidth',2);
        %plot(lambda*res.L500.y,res.L500.hP,'k','linewidth',2);    
        
                 
           %Limiting behaviour        
        
        
        plot(yP,1 + delta*(log(yP/lambda)+Cin),'r-.','linewidth',1.);     hold on;
        markI = ((yP*lambda<1e-4) & (yP*lambda>7e-6));        
        yPI = 5e-6*(1:100)';
        plot(yPI,1 + delta*(log(yPI/lambda)+Cin)...
                  - delta^2*((log(yPI/lambda)).^2+2*log(yPI/lambda)),'r:','linewidth',1.);     hold on;
                
        set(gca,'XScale','log');    
        
        ylim([0.9 2.5]);
        xlim([1e-8 1e2]);            

        xlabel('$x$','Interpreter','Latex');
        ylabel('$dh/dx$','Interpreter','Latex');
        
        set(gca,'XTick',[1e-8,1e-6,1e-4,1e-2,1e0,1e2]);        
        
        SaveFigure('ForcedWettingLubrication');        
    end