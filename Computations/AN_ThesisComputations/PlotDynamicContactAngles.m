function PlotDynamicContactAngles
    t = 0:0.01:pi;

    if(exist('D:\','dir'))
        dir = 'D://SyncGit/Projects/Data/ExperimentalContactAngle/';
    elseif(exist('/Users/NoldAndreas/','dir'))
        dir = '/Users/NoldAndreas/Documents/SyncGit/Projects/Data/ExperimentalContactAngle/';
    end
	close all;
    
    nData = 1;
    data  = {};
    
    L_Lambda    = [10^4,exp(1.92/3)];
    lines = {'k--','k:','k-.'};
    
    
    HoffmannFit = LoadHoffmannData('Fit','-',[]);
    
    nData = 1;

    data{nData} = LoadDataFluidData('Stroem/Stroem_ParaffinOil_untreatedPS','s','r',struct('rad','deg','Ca','V'));  
    %LoadDataFluidData('Stroem/Stroem_ParaffinOil_PTFE','d','r'); %(no zero static contact angle)
    data{nData} = LoadDataFluidData('Stroem/Stroem_SiliconeOil_I_untreatedPS','d','r',struct('rad','deg','Ca','V'));
    data{nData} = LoadDataFluidData('Stroem/Stroem_SiliconeOil_II_untreatedPS','o','r',struct('rad','deg','Ca','V'));
    data{nData} = LoadDataFluidData('Stroem/Stroem_SiliconeOil_II_oxidizedPS','v','r',struct('rad','deg','Ca','V'));
    
    data{nData} = LoadHoffmannData('BlackCircles','s','k');
    data{nData} = LoadHoffmannData('Triangles','d','k');
    data{nData} = LoadHoffmannData('Crosses','^','k');
    data{nData} = LoadHoffmannData('Hexagons','o','k');
    data{nData} = LoadHoffmannData('Squares','v','k');
    
    LoadDataFluidData('Kavehpour/eta_0_007','s','m',struct('rad','rad','Ca','Ca'));
    LoadDataFluidData('Kavehpour/eta_0_5','d','m',struct('rad','rad','Ca','Ca'));
    LoadDataFluidData('Kavehpour/eta_1_0','^','m',struct('rad','rad','Ca','Ca'));
    LoadDataFluidData('Kavehpour/eta_5_0','o','m',struct('rad','rad','Ca','Ca'));
    LoadDataFluidData('Kavehpour/eta_10_0','v','m',struct('rad','rad','Ca','Ca'));

	LoadDataFluidData('RameGaroff/CH3','v','g',struct('rad','deg','Ca','Ca'));
    LoadDataFluidData('RameGaroff/OH','o','g',struct('rad','deg','Ca','Ca'));
   
    PlotGThetaOverCa();
    PlotThetaOverCaPlusF(data,'FHoffmann');
    PlotThetaOverCaPlusF(data,'FGR');       
    
    PlotThetaOverCa(data);  
             
    
    function PlotThetaOverCa(data)
        figure('color','white','Position',[0 0 800 800]);
        
        for i = 1:length(data)             
            plotLegend{i} = data{i}.legend;
            if(isfield(data{i},'lw'))                
                semilogx(data{i}.Ca,data{i}.theta,...
                    data{i}.symbol,'linewidth',data{i}.lw); 
            else                
                semilogx(data{i}.Ca,data{i}.theta,...
                    data{i}.symbol,'MarkerFaceColor',data{i}.color,...
                    'MarkerSize',data{i}.MarkerSize); hold on;
            end
            hold on;
        end
        noPlots = i + 1;
        
        thetaM = (0:0.01:pi)';
        G_thM  = GHR_lambdaEta0(thetaM);

        for k = 1:length(L_Lambda)
            plot(G_thM/log(L_Lambda(k)),thetaM*180/pi,lines{k},'linewidth',2);
            plotLegend{noPlots} = ['analytical prediction with L/lambda = ',num2str(L_Lambda(k))];
            noPlots = noPlots + 1;
        end
        
        h_legend = legend(plotLegend,'Location','Northwest');          
        set(h_legend,'FontSize',10);        
        
       % xlim([1e-5 100]);        ylim([0 180]);
        xlabel('$Ca$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta_m [^\circ]$','Interpreter','Latex','fontsize',20);

        set(gca,'fontsize',20);

        print2eps('ThetaOverCa',gcf);  
        saveas(gcf,'ThetaOverCa.fig');
         
    end
    function PlotThetaOverCaPlusF(data,FHoffmann_FGHR)
        figure('color','white','Position',[0 0 800 800]);
        noPlots = 1;                        
        
        
        if(strcmp(FHoffmann_FGHR,'FHoffmann'))
            semilogx(HoffmannFit.Ca+HoffmannFit.Feq,HoffmannFit.theta,...
                            HoffmannFit.symbol,'linewidth',HoffmannFit.lw);         
            plotLegend{noPlots} = HoffmannFit.legend;
            noPlots = noPlots + 1;        
        elseif(strcmp(FHoffmann_FGHR,'FGR'))
            thetaM = (0:0.01:pi)';
            G_thM  = GHR_lambdaEta0(thetaM);    
            for k = 1:length(L_Lambda)
                semilogx(G_thM/log(L_Lambda(k)),thetaM*180/pi,lines{k},'linewidth',2);  hold on;
                plotLegend{noPlots} = ['Fit to G with L/lambda = ',num2str(L_Lambda(k))];
                noPlots = noPlots + 1;  
            end
        end  
        
        hold on;
                
        for i = 1:length(data)
            
            if(data{i}.thetaEq == 0)
                data{i}.Feq = 0;
            else
                if(strcmp(FHoffmann_FGHR,'FHoffmann'))                
                    feq = interp1q(HoffmannFit.theta,HoffmannFit.Ca,data{i}.thetaEq);                    
                elseif(strcmp(FHoffmann_FGHR,'FGR'))
                    feq = interp1q(thetaM*180/pi,G_thM/log(10^4),data{i}.thetaEq);                    
                end
                if(isfield(data{i},'Feq'))
                    disp(['Feq differs from data in file ',data{i}.legend,' :',num2str(abs(feq-data{i}.Feq))]);
                end
                data{i}.Feq = feq;
            end
            
            plotLegend{noPlots} = data{i}.legend;
            noPlots = noPlots + 1;
            if(isfield(data{i},'lw'))                
                semilogx(data{i}.Ca+data{i}.Feq,data{i}.theta,...
                    data{i}.symbol,'linewidth',data{i}.lw); 
            else                
                semilogx(data{i}.Ca+data{i}.Feq,data{i}.theta,...
                    data{i}.symbol,...
                    'color',data{i}.color,'MarkerFaceColor',data{i}.color,...
                    'MarkerSize',data{i}.MarkerSize); hold on;
            end                        
            hold on;
        end                      
        
        h_legend = legend(plotLegend,'Location','Northwest');          
        set(h_legend,'FontSize',10);        
        
        for i = 1:length(data)
            semilogx(data{i}.Feq,data{i}.thetaEq,...
                    data{i}.symbol,...
                    'LineWidth',2,...
                    'MarkerEdgeColor','g','MarkerFaceColor',data{i}.color,...
                    'MarkerSize',data{i}.MarkerSize+2); hold on; 
        end

        %xlim([1e-5 100]); ylim([0 180]);
        
        ylabel('$\theta_m [^\circ]$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',20);
        
         if(strcmp(FHoffmann_FGHR,'FHoffmann'))                
            xlabel('$Ca + F_H(\theta_{eq})$','Interpreter','Latex','fontsize',20);
            filename = 'ThetaOverCaF_H';
        elseif(strcmp(FHoffmann_FGHR,'FGR'))
            xlabel('$Ca + F_G(\theta_{eq})$','Interpreter','Latex','fontsize',20);
            filename = 'ThetaOverCaF_G';
         end        

        print2eps(filename,gcf);  
        saveas(gcf,[filename,'.fig']);
        
    end
    function PlotGThetaOverCa()
        
        figure('color','white','Position',[0 0 800 800]);
        noPlots = 1;                                
        
        for i = 1:length(data)
                        
            g   = GHR_lambdaEta0(data{i}.theta*pi/180) - GHR_lambdaEta0(data{i}.thetaEq*pi/180);
            
            plotLegend{noPlots} = data{i}.legend;
            noPlots = noPlots + 1;
            
            if(isfield(data{i},'lw'))                
                loglog(data{i}.Ca,g,...
                    data{i}.symbol,'linewidth',data{i}.lw); 
            else                
                loglog(data{i}.Ca,g,...
                    data{i}.symbol,'MarkerFaceColor',data{i}.color,...
                    'MarkerSize',data{i}.MarkerSize); hold on;
            end
            hold on;
        end
                      
        xP = 10.^(-6:0.1:1);%(0:1e-6:10);
        for k = 1:length(L_Lambda)
            plot(xP,xP*log(L_Lambda(k)),lines{k},'linewidth',2);
            plotLegend{noPlots} = ['analytical prediction with L/lambda = ',num2str(L_Lambda(k))];
            noPlots = noPlots + 1;
        end
        
        h_legend = legend(plotLegend,'Location','Northwest');          
        set(h_legend,'FontSize',10);                       
                 
        xlabel('$Ca$','Interpreter','Latex','fontsize',20);        
        ylabel('$G(\theta_m) - G(\theta_{eq})$','Interpreter','Latex','fontsize',20);

        set(gca,'fontsize',15);
        
        %xlim([1e-5,1e2]);
        %ylim([1e-4,1e2]);

        print2eps('GThetaOverCa',gcf);  
        saveas(gcf,'GThetaOverCa.fig');
        
    end        
   
    function dataL =  LoadHoffmannData(name,symbol,color)
        fid = fopen([dir,'Hoffmann/Hoffmann_',name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',4); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        dataL.legend = char(y{1});
        dataL.Feq = x{1}(1);
        dataL.thetaEq = x{1}(2);
        dataL.Ca  = x{1}(3:end) - dataL.Feq;
        dataL.theta = x{2}(3:end);                
        
        dataL.color  = color;
        dataL.symbol = symbol;
        
        if(isempty(color))
            dataL.lw    = 2;
            %semilogx(data.Ca,data.theta,symbol,'linewidth',2); hold on;
        else
            dataL.MarkerSize = 10;
            %semilogx(data.Ca,data.theta,symbol,'MarkerFaceColor',color,'MarkerSize',10); hold on;
        end
        
        data{nData} = dataL;
        nData       = nData + 1;                
    end
    function dataL = LoadDataFluidData(name,symbol,color,opts)
        fid = fopen([dir,name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',4); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        
        dataL.gamma = x{1}(1);
        dataL.eta = x{1}(2);
        dataL.rho = x{1}(3);
        
        if(strcmp(opts.rad,'rad'))
            dataL.theta = x{1}(4:end)*180/pi;            
        elseif(strcmp(opts.rad,'deg'))
            dataL.theta = x{1}(4:end);
        end
        if((nargin > 3) && strcmp(opts.Ca,'Ca'))
            dataL.Ca    = x{2}(4:end);
        elseif((nargin > 3) && strcmp(opts.Ca,'V'))
            dataL.Ca    = x{2}(4:end)*dataL.eta/dataL.gamma*10^-3;
        end
        
        dataL.thetaEq = dataL.theta(dataL.Ca == 0);        
        dataL.legend  = char(y{1});
        
        dataL.color  = color;
        dataL.symbol = symbol;
        
        if(isempty(color))
            dataL.lw = 2;            
        else
            dataL.MarkerSize = 10;            
        end
        
        data{nData} = dataL;
        nData       = nData + 1;
               
    end
%     function z = GHR(t)
%         z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
%             + 1i/2*(  dilog(1+exp(1i*t)) + ...
%                       dilog(exp(1i*t)) ) - sin(t)/2;
%         disp(['Max imaginary part: ',num2str(max(imag(z)))]);
%         z = real(z);
%     end
end