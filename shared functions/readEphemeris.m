function [kep] = readEphemeris(filename)
%
% PROTOTYPE:
%  [kep] = readEphemeris(filename);
% 
% DESCRIPTION:
%   The function creates a matrix of orbital parameters from a txt file
%       
% INPUT:
%   filename[-]     File's name containing the orbital parameters    [char]
%                                                                        
% OUTPUT:
%   kep[?x6]        Matrix containing the orbital parameters extracted from
%                   the file, with the following column order:
%                       - Semi-major axis [km]
%                       - Eccentricity [-]
%                       - Inclination [deg]
%                       - RAAN [deg]
%                       - Argument of pericentre [deg]
%                       - True anomaly [deg]
% 
% FUNCTIONS CALLED:
%   (none)
%
% AUTHORS:
%   Giuseppe Brentino, Nicolò Galletta, Roberto Pistone Nascone, Virginia
%   Di Biagio Missaglia, 2022
%--------------------------------------------------------------------------

ephemeris = readmatrix(filename); %extract data from file

%% find end of useful data

found = 0;
row = 3;

while ~found
    if isnan(ephemeris(row, 1)) 
        found = 1;
    end
    row = row+1;
end

row = row-2;

%% generate output matrix

kep=[ephemeris(3:row,12), ephemeris(3:row,3), ...
    ephemeris(3:row,5:7), ephemeris(3:row,11)];
