function ThesisNanoscale_Fig2_NumParaemeterStudy()  
    
    AddPaths('ThesisNanoscale');

	CLT = [];
    DoParameterStudy(0.453,[135 134]);
    DoParameterStudy(1.071,[60 61]);

    
    function DoParameterStudy(epw,alpha_deg)                 

        %Setup result file for this Job        
        filename   = ['ThesisNanoscale_Fig2_NumParaemeterStudy',getTimeStr(),'.txt'];
        Struct2File(filename,v2struct(epw,alpha_deg),['Computed at ',datestr(now)]);    

        close all;  
        
        config = ThesisNanoscale_GetStandardConfig([],[]);
        opts   = v2struct(alpha_deg,epw,config);
        st     = DataStorage('NumericalStudy',@NumericalStudyCL,opts,[],[]);        
        
        cols = {'k','b','m'};
        %Load configuration \Config_2014_8_18_11_32_15.686.mat    
        f1 = figure('Color','white','Position',[0 0 350 300]);    

        for i = 1:length(alpha_deg)
            plot(st.results15.y2{i},st.results15.theta{i},[cols{i},':']); hold on; 
            plot(st.results20.y2{i},st.results20.theta{i},[cols{i},':']); hold on; 
            plot(st.results25.y2{i},st.results25.theta{i},[cols{i},':']); hold on; 
            plot(st.results30.y2{i},st.results30.theta{i},[cols{i},':']); hold on; 
            %plot(y2_35,theta_35,'m','linewidth',1.5); hold on; 
            
            plot([5 40],alpha_deg(i)*[1 1],[cols{i},'--'],'linewidth',2);
        end
        
        xlim([5 37]);
        plot([5 40],CLT.alpha_YCA*180/pi*[1 1],'r--','linewidth',2);
        
        xlabel('$y_2$','Interpreter','Latex');
        ylabel('$\theta [^\circ]$','Interpreter','Latex');        

        f2 = figure('Color','white');
        CLT.PlotContourResults({});

        inset2(f1,f2,0.5,[0.45,0.55]);
        close(f2);  

        ChangeDirData();
        SaveFigure(['Fig2_epw_',num2str(epw)]);
    end

    function st = NumericalStudyCL(st,misc)        
        alpha_deg = st.alpha_deg;
        epw       = st.epw;
        
        PlotAreaCart       = struct('y1Min',-15,'y1Max',10,...
                                    'y2Min',0.5,'y2Max',18.5,...
                                    'zMax',4,...
                                    'N1',100,'N2',100); 
          
        for i = 1:length(alpha_deg)
            
            clear 'CLT'
            
            config = ThesisNanoscale_GetStandardConfig(alpha_deg(i),epw);
            config.optsNum.PlotAreaCart = PlotAreaCart;                                    
                        
            CLT = ContactLineHS(config);
            CLT.Preprocess();            

            [y2_15{i},theta_15{i}] = GetY2Theta(15);
            [y2_20{i},theta_20{i}] = GetY2Theta(20);
            [y2_25{i},theta_25{i}] = GetY2Theta(25);
            [y2_30{i},theta_30{i}] = GetY2Theta(30);
            %[y2_35{i},theta_35{i}] = GetY2Theta(35);
          end
        
          st.results15 = struct('y2',y2_15,'theta',theta_15);
          st.results20 = struct('y2',y2_20,'theta',theta_20);
          st.results25 = struct('y2',y2_25,'theta',theta_25);
          st.results30 = struct('y2',y2_30,'theta',theta_30);
    end
          
	function [y2,theta] = GetY2Theta(y2Max)
        CLT.optsNum.maxComp_y2 = y2Max;
        CLT.ComputeEquilibrium(struct('solver','Picard'));     
        %CLT.ComputeEquilibrium();  
        [y2,theta] = CLT.PlotInterfaceAnalysisY2([5 (y2Max+3)]);        
    end
end