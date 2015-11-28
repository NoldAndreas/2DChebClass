function ThesisNanoscale_DifferentShapes

    AddPaths('ThesisNumerics');
    
	%SavePic('SimulationAnnulus');    
    %SavePic('SimulationDisk');
    %SavePic('SimulationInfAnnulus');       
    %SavePic('SimulationSphere');    
    %SavePic('SimulationBall');        
    
    SavePic('SimulationStripMinusDisk');
    SavePic('SimulationHalfSpaceMinusHalfDisk');
    SavePic('SimulationHalfStripMinusDisk');    
    SavePic('SimulationWedge');
    SavePic('SimulationWedgeCut');
    SavePic('SimulationWedgeCutSide');
    

    function SavePic(strf)
        f = str2func(strf);
        
        [~,res] = f();
        fh = res.fig_handles{1};        
        set(0,'currentfigure', fh);        
        set(fh,'Position',[0 0 150 100],'color','white');                                                            
        ax = get(fh,'children');
                
        xlim = get(ax,'xlim');
        ylim = get(ax,'ylim');

        r = (ylim(2)-ylim(1))/(xlim(2)-xlim(1));
        pbaspect(ax,[(xlim(2)-xlim(1)) (ylim(2)-ylim(1)) 1]);
        
        %pbaspect(ax,[1 1 1]);
        xlabel(''); ylabel('');
        set(0, 'currentfigure',fh);
        SaveFigure(strf);
    end
end