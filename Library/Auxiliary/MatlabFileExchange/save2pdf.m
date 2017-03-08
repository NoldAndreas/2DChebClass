function save2pdf(pdfFileName,handle,dpi,saveEPS)

if(nargin<3)
    dpi=[];
end

if nargin<4
    saveEPS=false;
end

% Backup previous settings
prePaperType = get(handle,'PaperType');
prePaperUnits = get(handle,'PaperUnits');
preUnits = get(handle,'Units');
prePaperPosition = get(handle,'PaperPosition');
prePaperSize = get(handle,'PaperSize');

% Make changing paper type possible
set(handle,'PaperType','<custom>');

% Set units to all be the same
set(handle,'PaperUnits','inches');
set(handle,'Units','inches');

% Set the page size and position to match the figure's dimensions
position = get(handle,'Position');
set(handle,'PaperPosition',[0,0,position(3:4)]);
set(handle,'PaperSize',position(3:4));

% Save figure as an eps file
epsFileName=[pdfFileName(1:end-3) 'eps'];
if(isempty(dpi))
    print(handle,'-depsc2','-noui','-painters',epsFileName)
else
    print(handle,'-depsc2','-noui',['-r' num2str(dpi)],epsFileName)
end

if verLessThan('matlab','8.4')
    % Open file and read it in
    fid = fopen(epsFileName, 'r');
    str = fread(fid);
    str = char(str');
    fclose(fid);

    % Find where the line types are defined
    id   = strfind(str, '% line types:');

    % Get the part of the file before this point
    beforeDefns = str(1:id-1);

    % find the first '/' which defines the styles
    [~, restOfFile] = strtok(str(id:end), '/');
    % ~ should be: % line types: solid, dotted, dashed, dotdash

    % find the first % which delimits the end of the definitions
    [~, restOfFile] = strtok(restOfFile, '%');
    % ~ contains the definitions to be replaced.  Should be:
    % /SO { [] 0 setdash } bdef
    % /DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef
    % /DA { [6 dpi2point mul] 0 setdash } bdef
    % /DD { [.5 dpi2point mul 4 dpi2point mul 6 dpi2point mul 4
    %  dpi2point mul] 0 setdash } bdef

    % Define the new line styles
    lineDefns   = sprintf('%% line types: solid, dotted, dashed, long dashed\n');
    solidLine   = sprintf('/SO { [] 0 setdash } bdef\n');
    dotLine     = sprintf('/DO { [3 dpi2point mul 3 dpi2point mul] 0 setdash } bdef\n');
    dashedLine  = sprintf('/DA { [6 dpi2point mul 6 dpi2point mul] 0 setdash } bdef\n');
    dashdotLine = sprintf('/DD { [12 dpi2point mul 6 dpi2point mul] 0 setdash } bdef\n');

    % Construct the new file with the new line style definitions
    newText     = [beforeDefns, lineDefns, solidLine, dotLine, dashedLine, dashdotLine, restOfFile];

    % Write file with new line definitions
    fid = fopen(epsFileName, 'w');
    fprintf(fid, '%s', newText);
    fclose(fid);
end

% Set ghostscript options to convert to pdf
GSopts = [' -q -dNOPAUSE -dBATCH -dEPSCrop -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="' pdfFileName '" -f "' epsFileName '"'];


% Get gs command depending on OS
switch computer
		case {'MAC','MACI','MACI64'}			
            gs= '/usr/local/bin/gs';
		case {'PCWIN','PCWIN64'}
            gs= 'gswin32c.exe';
        otherwise
            gs= 'gs';
end

% Use gs to convert to pdf
eps2pdfCmd=[ gs GSopts];
system(eps2pdfCmd);

% Delete the eps file if it's not required
if(~saveEPS)
    delete(epsFileName);
end

% Restore the previous settings
set(handle,'PaperType',prePaperType);
set(handle,'PaperUnits',prePaperUnits);
set(handle,'Units',preUnits);
set(handle,'PaperPosition',prePaperPosition);
set(handle,'PaperSize',prePaperSize);
