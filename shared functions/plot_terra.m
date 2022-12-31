function plot_terra
%
% PROTOTYPE:
%   plot_terra;
%
% DESCRIPTION:
%  Creates a 3D model of the Earth loading an image for texture map from a
%  local directory.
%  
% INPUT:
%   (none)
%
% OUTPUT:
%   (none)
%
% FUNCTIONS CALLED:
%   (none)
%
% AUTHORS: 
%   Ryan Gray
%
% VERSIONS:
%   8 Sep 2004, 31 Jan 2006, 16 Oct 2013
%   2022 changed background and Earth image by: Nicolò Galletta, Virginia Di
%   Biagio Missaglia, Roberto Pistone Nascone, Giuseppe Brentino  
%--------------------------------------------------------------------------

%% Options
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
%GMST0 = []; % Don't set up rotatable globe (ECEF)
GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe image
image_file = 'Earth.png';

% Mean spherical earth
erad    = 6371.0087714; % equatorial radius (meters)
prad    = 6371.0087714; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec) 

%% Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

%% Texturemap the globe

% Load Earth image for texture map
cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
hold on;