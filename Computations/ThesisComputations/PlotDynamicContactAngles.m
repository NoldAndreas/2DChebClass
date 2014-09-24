function PlotDynamicContactAngles   

    if(exist('D:\','dir'))
        dir = 'D://SyncGit/Projects/Data/ExperimentalContactAngle/';
    elseif(exist('/Users/NoldAndreas/','dir'))
        dir = '/Users/NoldAndreas/Documents/SyncGit/Projects/Data/ExperimentalContactAngle/';
    end
	close all;
    
    nData = 1;
    
    HoffmannFit = LoadHoffmannData('Fit','-',[]);

%     data{nData} = LoadDataFluidData('Stroem/Stroem_ParaffinOil_untreatedPS','s','r');  
%     %LoadDataFluidData('Stroem/Stroem_ParaffinOil_PTFE','d','r'); %(no zero static contact angle)
%     data{nData} = LoadDataFluidData('Stroem/Stroem_SiliconeOil_I_untreatedPS','d','r');    
%     data{nData} = LoadDataFluidData('Stroem/Stroem_SiliconeOil_II_untreatedPS','o','r');    
%     data{nData} = LoadDataFluidData('Stroem/Stroem_SiliconeOil_II_oxidizedPS','v','r');    
%     
%     data{nData} = LoadHoffmannData('BlackCircles','s','k');
%     data{nData} = LoadHoffmannData('Triangles','d','k');
%     data{nData} = LoadHoffmannData('Crosses','^','k');
%     data{nData} = LoadHoffmannData('Hexagons','o','k');
%     data{nData} = LoadHoffmannData('Squares','v','k');            
    
    data{nData} = LoadDataFluidData('Kavehpour/eta_0_5','v','k'); 
    
   
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
        
        h_legend = legend(plotLegend,'Location','Northwest');          
        set(h_legend,'FontSize',10);        

        thetaM = (0:0.01:pi)';
        G_thM  = GHR(thetaM);

        plot(G_thM/log(10^4),thetaM*180/pi,'k--','linewidth',2);

        xlim([1e-5 100]);
        ylim([0 180]);
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
            G_thM  = GHR(thetaM);    
            semilogx(G_thM/log(10^4),thetaM*180/pi,'k--','linewidth',2); 
            plotLegend{noPlots} = 'Fit to G';
            noPlots = noPlots + 1;  
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

        xlim([1e-5 100]);
        ylim([0 180]);
        
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
                        
            g   = GHR(data{i}.theta*pi/180) - GHR(data{i}.thetaEq*pi/180);
            
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
                      
        xP = (0:1e-5:10);
        plot(xP,xP*log(1e4),'linewidth',2);
        plotLegend{noPlots} = 'analytical prediction with L/lambda = 10^4';
        noPlots = noPlots + 1;
        
        h_legend = legend(plotLegend,'Location','Northwest');          
        set(h_legend,'FontSize',10);                       
                 
        xlabel('$Ca$','Interpreter','Latex','fontsize',20);        
        ylabel('$G(\theta_m) - G(\theta_{eq})$','Interpreter','Latex','fontsize',20);

        set(gca,'fontsize',15);
        
        xlim([1e-5,1e2]);
        ylim([1e-4,1e2]);

        print2eps('GThetaOverCa',gcf);  
        saveas(gcf,'GThetaOverCa.fig');
        
    end        
   
    function data =  LoadHoffmannData(name,symbol,color)
        fid = fopen([dir,'Hoffmann/Hoffmann_',name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',4); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        data.legend = char(y{1});
        data.Feq = x{1}(1);
        data.thetaEq = x{1}(2);
        data.Ca  = x{1}(3:end) - data.Feq;
        data.theta = x{2}(3:end);                
        
        data.color  = color;
        data.symbol = symbol;
        
        if(isempty(color))
            data.lw    = 2;
            %semilogx(data.Ca,data.theta,symbol,'linewidth',2); hold on;
        else
            data.MarkerSize = 10;
            %semilogx(data.Ca,data.theta,symbol,'MarkerFaceColor',color,'MarkerSize',10); hold on;
        end
        
        nData = nData + 1;
                
    end
    function data =  LoadDataFluidData(name,symbol,color)
        fid = fopen([dir,name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',4); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        
        data.gamma = x{1}(1);
        data.eta = x{1}(2);
        data.rho = x{1}(3);
        
        data.theta = x{1}(4:end);
        data.Ca    = x{2}(4:end)*data.eta/data.gamma*10^-3;
        
        data.thetaEq = data.theta(data.Ca == 0);        
        data.legend  = char(y{1});
        
        data.color  = color;
        data.symbol = symbol;
        
        if(isempty(color))
            data.lw = 2;            
        else
            data.MarkerSize = 10;            
        end
        
        nData = nData + 1;
               
    end
    function z = GHR(t)
        z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
            + 1i/2*(  dilog(1+exp(1i*t)) + ...
                      dilog(exp(1i*t)) ) - sin(t)/2;
        disp(['Max imaginary part: ',num2str(max(imag(z)))]);
        z = real(z);
    end
end