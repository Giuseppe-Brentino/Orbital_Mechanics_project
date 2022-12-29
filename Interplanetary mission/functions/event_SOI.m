function [value, isterminal, direction] = event_SOI(t, s, settings, SOI, varargin)
%
% DESCRIPTION:
%   This function checks at every step of the ODE if the s/c is into the
%   SOI of Saturn. It ends the integration as soon as the s/c exits from
%   Saturn's SOI.
%
%  INPUT :
%   t[?,1]      time vector from ODE                    [s]
%   s[?,6]      state matrix from ODE contains          [km; km/s]
%   settings    struct containing the following parameters: 
%                   - settings.mu[1]            planetary constant [km^3/s^2]
%                   - settings.perturbations    logical value, activates
%                                               perturbations
%   varargin    cell array containig optional input arguments          
%
%  OUTPUT:
%	value[1]    expression that checks if the s/c is inside the SOI of the
%	            planet. If its numerical value is equal to 0 it activates
%	            "isterminal".                           [km]
%   isterminal  flag that terminates the integration when value=0
%   direction   default flag that keeps the integration going
%
%  FUNCTIONS CALLED:
%   (none)
%
% AUTHORS:
%   Virginia Di Biagio Missaglia, Roberto Pistone Nascone, Giuseppe
%   Brentino, Nicolò Galletta 2022
%
% -------------------------------------------------------------------------

value = norm(s(1:3)) - SOI;

isterminal = 1;

direction = 0;

end