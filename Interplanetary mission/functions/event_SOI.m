function [value, isterminal, direction] = event_SOI(t, s, settings,SOI, varargin)

value = norm(s(1:3)) - SOI;

isterminal = 1;

direction = 0;

end