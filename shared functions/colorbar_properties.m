function obj = colorbar_properties(obj,varargin)
% Change the properties of the colorbar named 'obj' to match the style of
% the figure.
%  
%  INPUT:
%  obj: object representing the colorbar
%  varargin: Use another variable to add a title to the colorbar
%  
% CONTRIBUTORS:
% Giuseppe Brentino, Roberto Pistone Nascone

fig_int=get(0,'defaultTextInterpreter');
obj.TickLabelInterpreter = fig_int;

fig_size = get(0,'DefaultAxesFontSize');
obj.FontSize = fig_size;

if nargin == 2
    name = varargin{1};
    title(obj,name,'Interpreter','latex');
end