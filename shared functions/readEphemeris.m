function [kep] = readEphemeris(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that creates a matrix of orbital parameters from a file
%       
% INPUT:
%   filename: variable of type char, it's the name of the file
%   containing the orbital parameters
%                                                                         
% OUTPUT:
%   kep:matrix containing the orbital parameters extracted from the file,
%   with the following order: 
%   [major semi-axis; eccentricity; inclination; RAAN; pericenter's anomaly; true anomaly]
% 
% Contributors:
%   Giuseppe Brentino, Nicolò Galletta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ephemeris=readmatrix(filename); % extract data from file

%% find end of useful data

found=0;
row=3;

while ~found
    if isnan(ephemeris(row,1)) 
        found=1;
    end
    row=row+1;
end

row=row-2;

%% generate output matrix

kep=[ephemeris(3:row,12), ephemeris(3:row,3), ...
    ephemeris(3:row,5:7), ephemeris(3:row,11)];
