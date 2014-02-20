function [N1,N2,PhysArea,SubArea,plotTimes] = LoadNumData(optsNum)

        PhysArea = optsNum.PhysArea;
        
        if(isfield(optsNum.PhysArea,'N1'))
            N1          = optsNum.PhysArea.N1; 
            N2          = optsNum.PhysArea.N2;
            PhysArea.N  = [N1;N2];
        else
            N1=optsNum.PhysArea.N(1); N2=optsNum.PhysArea.N(2);
            PhysArea.N1 = PhysArea.N(1);
            PhysArea.N2 = PhysArea.N(2);
        end
                      
        if(isfield(optsNum,'SubArea'))
            SubArea   = optsNum.SubArea;           
            SubArea.N = [SubArea.N1;SubArea.N2];
        else 
            SubArea   = PhysArea;
        end

        if(nargout == 5)
            plotTimes  = optsNum.plotTimes;        
        end
        
        %SubArea Limits
        if((SubArea.y1Max < PhysArea.y1Min) || (SubArea.y1Max > PhysArea.y1Max))
             err = MException('LoadNumData:SubAreautOfRange', ...
                    'SubArea.y1Max value is outside expected range');
             throw(err);
        end
        if((SubArea.y1Min < PhysArea.y1Min) || (SubArea.y1Min > PhysArea.y1Max))
             err = MException('LoadNumData:SubAreautOfRange', ...
                    'SubArea.y1Min value is outside expected range');
             throw(err);
        end
        if((SubArea.y2Max < PhysArea.y2Min) || (SubArea.y2Max > PhysArea.y2Max))
             err = MException('LoadNumData:SubAreautOfRange', ...
                    'SubArea.y2Max value is outside expected range');
             throw(err);
        end
        if((SubArea.y2Min < PhysArea.y2Min) || (SubArea.y2Min > PhysArea.y2Max))
             err = MException('LoadNumData:SubAreautOfRange', ...
                    'SubArea.y2Min value is outside expected range');
             throw(err);
        end        
        
        %Plot Ara Limits
        if((optsNum.PlotArea.y1Min < PhysArea.y1Min) || (optsNum.PlotArea.y1Min > PhysArea.y1Max))
             err = MException('LoadNumData:PlotAreaOutOfRange', ...
                    'PlotArea.y1Min value is outside expected range');
             throw(err);
        end
        if((optsNum.PlotArea.y1Max < PhysArea.y1Min) || (optsNum.PlotArea.y1Max > PhysArea.y1Max))
             err = MException('LoadNumData:PlotAreaOutOfRange', ...
                    'PlotArea.y1Max value is outside expected range');
             throw(err);
        end
        if((optsNum.PlotArea.y2Min < PhysArea.y2Min) || (optsNum.PlotArea.y2Min > PhysArea.y2Max))
             err = MException('LoadNumData:PlotAreaOutOfRange', ...
                    'PlotArea.y2Min value is outside expected range');
             throw(err);
        end
        if((optsNum.PlotArea.y2Max < PhysArea.y2Min) || (optsNum.PlotArea.y2Max > PhysArea.y2Max))
             err = MException('LoadNumData:PlotAreaOutOfRange', ...
                    'PlotArea.y2Max value is outside expected range');
             throw(err);
        end        
    
end