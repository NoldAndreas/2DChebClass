function optsLegOut=makeLegLine(optsLeg)
%optsLegOut=makeLegLine(optsLeg)
%  Returns a structure of legend and line data, such as legend labels, line
%  styles and colours.  In particular, it strips out data for the
%  simulations we're not running (false in doStoc/doDDFT) and only includes
%  type 'rv' simulations in momentum plots.  Also adds equilibrium legends.
% INPUTS:
%  optsLeg -- a structure containg
%         doStoc       (whether to do each of the stochastic simulations, 
%                        of the form {true, false,...})
%         doDDFT       (same for DDFT)
%         sampleFinal  (true/false whether we sample, and hence plot, 
%                        the final eq)
%         stocName     (legend labels for each stochastic simulation, 
%                        of the form {name1,name2,...})
%         DDFTName     (same for DDFT)
%         stocStyle    (line styles for plotting stochastic data, 
%                        of the same form)
%         DDFTStyle    (same for DDFT)
%         stocColour   (line colours for plotting stochastic data, 
%                        of the same form)
%         DDFTColour   (same for DDFT)
%         stocType     (type of stochastic simulations, {'rv','r',...})
%         DDFTType     (same for DDFT)
%         legPos       (legend position, e.g. 'NorthWest')
%         oneLeg       (either 'off' or anything else, to plot a single 
%                        legend at the top of the figure)
%         perRow       (integer max no. labels per row in single legend)
%
% OUTPUTS:
%  optsLegOut -- a structure of size 1x3 containing, for General, Initial
%                and Final plots
%             legTextR          (legend text for position plots,
%                                 as stocName)
%             legTextP          (same for momentum plots, this has
%                                 no type 'r calculations)
%             lineStyleStoc     (stochastic line styles, as input)
%             lineStyleDDFT     (same for DDFT) 
%             lineColourStoc    (stochastic line colours, as input)
%             lineColourDDFT    (same for DDFT)
%             stocType          (as input)
%             DDFTType          (as input)
%             lineStyle         (=[] ; dummy for use later)
%             lineColour        (=[] ; dummy for use later)
%             legText           (=[] ; dummy for use later)
%             legPos            (as input)
%             oneLeg            (as input)
%             perRow            (as input)


% get input options
% doStoc=optsLeg.doStoc;
% anyStoc=any(cell2mat(doStoc));
% doDDFT=optsLeg.doDDFT;
% anyDDFT=any(cell2mat(doDDFT));

anyStoc=optsLeg.anyStoc;
anyDDFT=optsLeg.anyDDFT;

if(anyStoc)
    doStoc=optsLeg.doStoc;
else
    doStoc=[];
end

if(anyDDFT)
    doDDFT=optsLeg.doDDFT;
else
    doDDFT=[];
end

% set default text for the legend
legTextR={};
legTextP=legTextR;

% get equilibrium styles
eqStyle=optsLeg.eqStyle;
eqColour=optsLeg.eqColour;
eqMarker=optsLeg.eqMarker;
eqText=optsLeg.eqText;

% set up stochastic legend information
if(anyStoc)
    
    % get input information
    sampleFinal=optsLeg.sampleFinal;
    %stocName=optsLeg.stocName;
    stocTitle=optsLeg.stocTitle;
    stocStyle=optsLeg.stocStyle;
    stocColour=optsLeg.stocColour;
    stocMarker=optsLeg.stocMarker;
    stocType=optsLeg.stocType;
    
    % set up placeholders
    lineStyleStoc={};
    lineMarkerStoc=lineStyleStoc;
    lineColourStoc=lineStyleStoc;

    % add in legend and line style for each of the calculations we're doing for
    % each of the General (1), Initial (2) and Final (3) cases
    for iStoc=1:length(doStoc)
        % only do for those stochastic simulations we want to run
        if(doStoc{iStoc})
            % add in legend text for position plots
            nSpecies=length(stocTitle{iStoc});
            % do for each species within the stochastic run
            for iSpecies=1:nSpecies            
                legTextR=cat(2,legTextR,stocTitle{iStoc}(iSpecies));
                if(strcmp(stocType{iStoc},'rv'))
                    % and for momentum plots if required
                    legTextP=cat(2,legTextP,stocTitle{iStoc}(iSpecies));
                end
            end
            
            % add in the line style and colour data
            lineStyleStoc=cat(2,lineStyleStoc,stocStyle(iStoc));
            lineMarkerStoc=cat(2,lineMarkerStoc,stocMarker(iStoc));
            lineColourStoc=cat(2,lineColourStoc,stocColour(iStoc));

        end
    end
  
end

if(anyDDFT)
    % same for DDFT    
    DDFTTitle=optsLeg.DDFTTitle;
    DDFTStyle=optsLeg.DDFTStyle;
    DDFTColour=optsLeg.DDFTColour;
    DDFTMarker=optsLeg.DDFTMarker;
    DDFTType=optsLeg.DDFTType;

    DDFTlegTextR={};
    DDFTlegTextP={};
    lineStyleDDFT={};
    lineMarkerDDFT=lineStyleDDFT;
    lineColourDDFT=lineStyleDDFT;

    for iDDFT=1:length(doDDFT)
        % add in legend text for each DDFT calculation we're doing
        if(doDDFT{iDDFT})
            
            nSpecies=length(DDFTTitle{iDDFT});
            
            for iSpecies=1:nSpecies
                % position text
                DDFTlegTextR=cat(2,DDFTlegTextR,DDFTTitle{iDDFT}(iSpecies));
                if(strcmp(DDFTType{iDDFT},'rv'))
                    % momentum text if needed
                    DDFTlegTextP=cat(2,DDFTlegTextP,DDFTTitle{iDDFT}(iSpecies));
                end
                
            end
            % style and colour
            lineStyleDDFT=cat(2,lineStyleDDFT,DDFTStyle(iDDFT));
            lineMarkerDDFT=cat(2,lineMarkerDDFT,DDFTMarker(iDDFT));
            lineColourDDFT=cat(2,lineColourDDFT,DDFTColour(iDDFT));
        end
    end
    
end


% strip out any stoc calculations we don't do, and add in equilibrium info
if(anyStoc)
    % stocTypePlot contains the types only for the calculations we want to do
    % (determined by doStoc)
    stocTypePlot=stocType(cell2mat(doStoc));
    
end

% strip out DDFT calculations we don't do
if(anyDDFT)
    DDFTTypePlot=DDFTType(cell2mat(doDDFT));
end
    
if(anyStoc && anyDDFT)

    optsLegOut.legTextR = cat(2,legTextR,DDFTlegTextR);
    optsLegOut.legTextP = cat(2,legTextP,DDFTlegTextP);
    
elseif(anyStoc)
    
    optsLegOut.legTextR = legTextR;
    optsLegOut.legTextP = legTextP;
    
else

    optsLegOut.legTextR = DDFTlegTextR;
    optsLegOut.legTextP = DDFTlegTextP;

end

if(anyStoc)
    optsLegOut.lineStyleStoc=lineStyleStoc;
    optsLegOut.lineMarkerStoc=lineMarkerStoc;
    optsLegOut.lineColourStoc=lineColourStoc;
    optsLegOut.stocType=stocTypePlot;
end

if(anyDDFT)
    optsLegOut.lineStyleDDFT=lineStyleDDFT;
    optsLegOut.lineMarkerDDFT=lineMarkerDDFT;
    optsLegOut.lineColourDDFT=lineColourDDFT;
    optsLegOut.DDFTType=DDFTTypePlot;    
end

optsLegOut.lineStyle=[];
optsLegOut.lineColour=[];
optsLegOut.legText=[];
optsLegOut.legPos=optsLeg.legPos;
optsLegOut.oneLeg=optsLeg.oneLeg;
optsLegOut.perRow=optsLeg.perRow;

optsLegOut.eqStyle=eqStyle;
optsLegOut.eqMarker=eqMarker;
optsLegOut.eqColour=eqColour;



