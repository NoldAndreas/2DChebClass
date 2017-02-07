function optsPlotV = getOptsPlotV(optsStruct)

temp=optsStruct.rMin;
dim=length(temp);

if(dim==1)
    YPlot=cell(1);   % Cartesian
    YPlotV=cell(1);  % Geometry
elseif(dim==2)
    YPlot1=cell(1);   % Cartesian
    YPlot2=cell(1);   % Cartesian
    YPlot1V=cell(1);  % Geometry
    YPlot2V=cell(1);  % Geometry
end

nSpecies=optsStruct.nSpecies;

temp=struct('rMin',optsStruct.rMin,'rMax',optsStruct.rMax);

rMin=temp.rMin;
rMax=temp.rMax;


if(dim==1)
    NPlot=optsStruct.NPlot-1;

    rRange=rMax-rMin;
    yPlot=rMin:rRange/NPlot:rMax;
    YPlot=yPlot(:);
    YPlotV=YPlot;

elseif(dim==2)
    NPlot1=optsStruct.NPlot1-1;
    NPlot2=optsStruct.NPlot2-1;

    rRange1=rMax{1}-rMin{1};
    rRange2=rMax{2}-rMin{2};
    yPlot1=rMin{1}:rRange1/NPlot1:rMax{1};
    yPlot2=rMin{2}:rRange2/NPlot2:rMax{2};
    yPlot1=yPlot1(:);
    yPlot2=yPlot2(:);

    [YPlot1,YPlot2]=meshgrid(yPlot1,yPlot2);

    switch optsStruct.geom

        case 'planar2D'
            YPlot1Vtemp=YPlot1;
            YPlot2Vtemp=YPlot2;

        case 'polar2D'
            [Theta,R]=cart2pol(YPlot1,YPlot2);
            YPlot1Vtemp=R(:);
            YPlot2Vtemp=Theta(:);
            %YPlot2Vtemp=mod(Theta,2*pi);
    end

    YPlot1V=repmat(YPlot1Vtemp(:),1,nSpecies);
    YPlot2V=repmat(YPlot2Vtemp(:),1,nSpecies);

end

if(dim==1)
    optsPlotV=struct('YPlot',YPlot,'YPlotV',YPlotV, ...
                     'NPlot',optsStruct.NPlot,'V1DV1',optsStruct.V1DV1);
elseif(dim==2)
    optsPlotV=struct('YPlot1',YPlot1,'YPlot2',YPlot2, ...
                      'YPlot1V',YPlot1V,'YPlot2V',YPlot2V, ...
                      'NPlot1',optsStruct.NPlot1,'NPlot2',optsStruct.NPlot2, ...
                      'quiverSkip',optsStruct.quiverSkip,'V1DV1',optsStruct.V1DV1);
end

end