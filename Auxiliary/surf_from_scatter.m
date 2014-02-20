function surf_from_scatter(x,y,z)
% %Copyright (c) 2009, The MathWorks, Inc.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

    %% Making Surface Plots From Scatter Data
    % How do you turn a collection of XYZ triplets into a surface plot? This is
    % the most frequently asked 3D plotting question that I got when I was in
    % Tech Support.

    %% Load the data

    %load seamount
    %who -file seamount

    %%
    % The problem is that the data is made up of individual (x,y,z)
    % measurements. It isn't laid out on a rectilinear grid, which is what the
    % SURF command expects. A simple plot command isn't very useful.

    %%plot3(x,y,z,'.-')

    %% Little triangles
    % The solution is to use Delaunay triangulation. Let's look at some
    % info about the "tri" variable.

    tri = delaunay(x,y);
    %plot(x,y,'.');

    %%
    % How many triangles are there?

    %[r,c] = size(tri);
    %disp(r)

    %% Plot it with TRISURF

    h = trisurf(tri, x, y, z);
    axis vis3d

    %% Clean it up

    %axis off
    %l = light('Position',[-50 -15 29])
    %set(gca,'CameraPosition',[208 -50 7687])
    lighting phong
    shading interp
    colorbar EastOutside
end