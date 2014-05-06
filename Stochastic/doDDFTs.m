function DDFTStructFull=doDDFTs(optsNumFull,optsPhysFull,optsPlot)
% DDFTStruct=doDDFTs(optsNumBoth,optsPhysGIF) 
%   determines which DDFT calculations are desired and performs them.
%
% INPUTS: 
%  optsNumFull -- a structure of size [1 nStoc] containing at least:
%         doDDFT          (true/false, whether to do DDFT calculation)
%         loadDDFT        (true/false, whether to load DDFT data)
%         saveFileDDFT    (string, save file name)
%         saveDDFT        (true/false, whether to save DDFT data)
%         plotTimes       (vector of times at which to save data)
%         type            (should be one of 'rv','r')
%         DDFTDir         (directory containing code)
%         DDFTCode        (name of DDFT .m file)
%         [DDFTOpts]    (simulation-specific options)
%
%  optsPhysGIF -- a structure of size [1 3] containing at least:
%         geom          (geometry; 'spherical'/'planar')
%         sigma         (particle diameter)
%         dim           (dimension; should always be 1)
%         D0            (diffusion constant)
%         gamma         (friction)
%         kBT           (temperature)
%         mu            (initial guess for mu)
%         nParticles    (number of particles)
%         V1DV1         (strings identifying potential and gradient functions
%                        {'General','Initial','Final'})
%         [pot params]  (parameters for potential where
%                           optsPhysGIF(1)=dynamic
%                           optsPhysGIF(2)=initial
%                           optsPhysGIF(3)=final)
%
% OUTPUTS:
%  DDFTStruct -- a structure of size [1 nDDFT] containing at least
%   rhov        -- a matrix of size (2 nPoints,nPoints) with colums of the
%                 form [rho; v] at times plotTimes.
%   r           -- a vector of size (nPoints, 1) giving Chebyshev points
%   w           -- a vector of size (1, nPoints) giving Chebyshev weights
%   meanR       -- mean positions over time
%   meanP       -- mean momenta over time

% number of DDFT calculations we're doing
nDDFT=length(optsNumFull);

lineStyles   = optsPlot.lineStyleDDFT;
lineMarkers  = optsPlot.lineMarkerDDFT;
lineColours  = optsPlot.lineColourDDFT;

% do each calculation
for iDDFT=1:nDDFT
    % get relevant options
    optsDDFT=optsNumFull(iDDFT);
    optsPhys=optsPhysFull(iDDFT);
    
    optsPlot.lineStyleDDFT=lineStyles{iDDFT};
    optsPlot.lineMarkerDDFT=lineMarkers{iDDFT};
    optsPlot.lineColourDDFT=lineColours{iDDFT};
    
    fprintf(1,['Starting DDFT ' num2str(iDDFT) '/' num2str(nDDFT) ...
                ': ' optsDDFT.DDFTName ' ... ']);

    loadDDFT = optsDDFT.loadDDFT;
    optsDDFT = rmfield(optsDDFT,'loadDDFT');
    optsDDFT = rmfield(optsDDFT,'saveDDFT');
	optsDDFT = rmfield(optsDDFT,'doDDFT');
    optsDDFT = rmfield(optsDDFT,'DDFTName');
    
    opts.optsPhys = optsPhys;
    opts.optsNum  = optsDDFT;
            
    DDFTdir = [optsPhys.potNames filesep 'DDFT' filesep optsDDFT.DDFTCode];            
    
    DDFTStruct = DataStorage(DDFTdir,@DDFTwrapper,opts,optsPlot,~loadDDFT);
                
    % save to full structure
    DDFTStructFull(iDDFT)=orderfields(DDFTStruct);  %#ok
    
%     % add in species coordinate for 1-species 2D DDFTs
%     if(optsPhysGIFFull(1).dim==2)
%         rho_t=DDFTStructFull(iDDFT).rho_t;
%         flux_t=DDFTStructFull(iDDFT).flux_t;
%         if(ndims(rho_t)==2)
%             temp=zeros(size(rho_t,1),1,size(rho_t,2));
%             temp(:,1,:)=rho_t;
%             DDFTStructFull(iDDFT).rho_t=temp;  %#ok
% 
%             temp=zeros(size(flux_t,1),1,size(flux_t,2));
%             temp(:,1,:)=flux_t;
%             DDFTStructFull(iDDFT).flux_t=temp; %#ok
%         end
%     end


    fprintf(1,'Finished\n');
end

    function output = DDFTwrapper(opts,optsPlot)
        optsP = opts.optsPhys;
        optsN  = opts.optsNum;
        f = str2func(optsN.DDFTCode);
        output = f(optsP,optsN,optsPlot);     
    end


end % doDDFTs