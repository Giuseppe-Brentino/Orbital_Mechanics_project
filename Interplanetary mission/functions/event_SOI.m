function [value, isterminal, direction] = event_SOI(t, s, settings, SOI, varargin)
%
% DESCRIPTION:
%   This function checks at every step of the ODE if the s/c is into the
%   SOI of Saturn. It ends the integration as soon as the s/c exits from
%   Saturn's SOI.
%
%  INPUT :
%   u[3x1]      Departure, flyby and arrival dates          [MJD2000]
%   i_dep[1]    Index of the departing planet to be used in the function
%               uplanet
%   i_fb[1]     Index of the flyby planet to be used in the function
%               uplanet
%   i_arr[1]    Index of the arriving planet to be used in the function
%               ephNEO
%
%  OUTPUT:
%	value[1]    Sum of three delta velocity: the one needed to perform the 
%               injection manoeuvre, the second one needed to perform the 
%               powered flyby and the third one to perform the arrival
%               manoeuvre.                                          [km/s]
%
%  FUNCTIONS CALLED:
%   (none)
%
% AUTHORS:
%   Virginia Di Biagio Missaglia, Roberto Pistone Nascone, Giuseppe
%   Brentino, Nicolò Galletta
%
% -------------------------------------------------------------------------

value = norm(s(1:3)) - SOI;

isterminal = 1;

direction = 0;

end