function [h_main, h_inset]=inset2(main_handle, inset_handle,inset_size,pos)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
% 
% An examle can found in the file: inset_example.m
% 
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.35;
end
if nargin < 4
    pos = [.7 .5];
end

inset_size=inset_size*.7;


inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,main_handle);
ax=get(main_handle,'Position');
if(length(inset_size) == 1)
    inset_size = inset_size*[1,1];
end
set(h_inset,'Position', [pos(1) pos(2) inset_size(1) inset_size(2)])