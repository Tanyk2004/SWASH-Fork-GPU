function plotunswash(matfile,basename,spar,type);
% Plots a SWASH output parameter on triangular mesh
%
% This parameter may either be defined per vertex (by default) or
% per cell
%
% Example:
%
% SWASH generated a binary Matlab file called 'a61dam03.mat'
% TRIANGLE generated files with basename 'a61dambr', e.g. 'a61dambr.ele'
% To make a plot of the water level, stored in each vertex, give the
% following command in Matlab:
%
%     plotunswash('a61dam03','a61dambr','Watlev','vertex')
%
%  or
%
%     plotunswash('a61dam03','a61dambr','Watlev')
%
% However if the water level is stored in each cell, then specify
% the following command line in Matlab:
%
%     plotunswash('a61dam03','a61dambr','Watlev','cell')
%
% For other output parameters, type the following command:
%
%    who -file a61dam03
%
%
% Author  : Marcel Zijlema
% Date    : March 01, 2022
% Version : 2.0

if nargin==3
   if ~exist('type','var')
      % type was not specified, so default is vertex
      type='node';
   end
elseif nargin~=4
   error('Wrong number of arguments. See "help plotunswash"')
end

eval(['load ' matfile]);                   % load binary file containing SWASH results
                                           % obtained using BLOCK command with COMPGRID-set
nodefile=[basename '.node'];
fid = fopen(nodefile);                     % load TRIANGLE node file
[nnode] = fscanf(fid,'%i',[1 4]);          % get number of nodes
ncol = 3+nnode(3)+nnode(4);                % specify number of columns in node file
nod = fscanf(fid,'%f',[ncol nnode(1)]);    % get nodes
x=nod(2,:); y=nod(3,:);                    % get coordinates of nodes
elefile=[basename '.ele'];
fid = fopen(elefile);                      % load TRIANGLE element file
[nelem] = fscanf(fid,'%i',[1 3]);          % get number of elements
ncol = 4+nelem(3);                         % specify number of columns in element file
tri = fscanf(fid,'%i',[ncol nelem(1)]);    % get elements
z=eval([spar]);                            % get SWASH output parameter
nz=size(z,2);                              % data size for checking
%
if strcmpi(type,'element') | strcmpi(type,'cell')
   if nz~=nelem(1)
      error('Output not consistent with mesh object')
   end
   t(1:3,:)=tri(2:4,:);                    % list of all triangles
   c=[x(t(:));y(t(:))];                    % list of coordinates of all triangles
   data=[z;z;z];                           % triplicate data so that each triangle has its own copy
   T=reshape(1:3*nelem(1),[3 nelem(1)]);   % connectivity table
   trisurf(T',c(1,:),c(2,:),data(:))       % make plot using trisurf
elseif strcmpi(type,'node') | strcmpi(type,'vertex')
   if nz~=nnode(1)
      error('Output not consistent with mesh object')
   end
   T=tri(2:4,:);                           % connectivity table
   trisurf(T',x,y,z)                       % make plot using trisurf
else
   error('Wrong type of mesh objects. See "help plotunswash"')
end
view(0,90);shading interp;                 % make 2D view and smooth plot
colormap(jet);colorbar;axis equal          % include colorbar and equal axes
