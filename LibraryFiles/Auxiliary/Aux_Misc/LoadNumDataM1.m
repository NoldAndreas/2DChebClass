function [N1,N2,PhysArea,SubArea,x1Plot,x2Plot,plotTimes] = LoadNumDataM1(optsNum,Maps)

        N1=optsNum.PhysArea.N1; N2=optsNum.PhysArea.N2;
        PhysArea = optsNum.PhysArea;
        
        if(~isfield(PhysArea,'L1') && (PhysArea.y1Max < inf))
            PhysArea.L1 = (PhysArea.y1Max - PhysArea.y1Min)/2;
        end

        if(~isfield(PhysArea,'L2') && (PhysArea.y2Max < inf))
            PhysArea.L2 = (PhysArea.y2Max - PhysArea.y2Min)/2;
        end
                      
        if(isfield(optsNum,'SubArea'))
            SubArea   = optsNum.SubArea;           
        else 
            SubArea   = PhysArea;
        end
        if(~isfield(SubArea,'L1') && (SubArea.y1Max < inf))
            SubArea.L1 = (SubArea.y1Max - SubArea.y1Min)/2;
        end

        if(~isfield(SubArea,'L2') && (SubArea.y2Max < inf))
            SubArea.L2 = (SubArea.y2Max - SubArea.y2Min)/2;
        end

        if(nargout == 7)
            plotTimes  = optsNum.plotTimes;        
        end
        
        CompSpace = Maps.CompSpace;
        
        x1Min = min(CompSpace(optsNum.PlotArea.y1Min,optsNum.PlotArea.y2Min,PhysArea),...
                    CompSpace(optsNum.PlotArea.y1Min,optsNum.PlotArea.y2Max,PhysArea));
        x1Max = max(CompSpace(optsNum.PlotArea.y1Max,optsNum.PlotArea.y2Min,PhysArea),...
                    CompSpace(optsNum.PlotArea.y1Max,optsNum.PlotArea.y2Max,PhysArea));                
        x2Min = min(CompSpace(optsNum.PlotArea.y1Min,optsNum.PlotArea.y2Min,PhysArea),...
                    CompSpace(optsNum.PlotArea.y1Max,optsNum.PlotArea.y2Min,PhysArea));
        x2Max = max(CompSpace(optsNum.PlotArea.y1Min,optsNum.PlotArea.y2Max,PhysArea),...
                    CompSpace(optsNum.PlotArea.y1Max,optsNum.PlotArea.y2Max,PhysArea));
        
        x1Plot = GetArray(x1Min,x1Max,optsNum.PlotArea.N1);
        x2Plot = GetArray(x2Min,x2Max,optsNum.PlotArea.N2);
        
        %Check correctness
        
        %PhysArea L
        if(PhysArea.L1 >= (PhysArea.y1Max - PhysArea.y1Min))
             err = MException('LoadNumData:LOutOfRange', ...
                    'PhysArea.L1 value is outside expected range');
             throw(err);
        end
        if(PhysArea.L2 >= (PhysArea.y2Max - PhysArea.y2Min))
             err = MException('LoadNumData:LOutOfRange', ...
                    'PhysArea.L2 value is outside expected range');
             throw(err);
        end
        
        %SubArea L
%         if(SubArea.L1 >= (SubArea.y1Max - SubArea.y1Min))
%              err = MException('LoadNumData:LOutOfRange', ...
%                     'SubArea.L1 value is outside expected range');
%              throw(err);
%         end
%         if(SubArea.L2 >= (SubArea.y2Max - SubArea.y2Min))
%              err = MException('LoadNumData:LOutOfRange', ...
%                     'SubArea.L2 value is outside expected range');
%              throw(err);
%         end
        
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