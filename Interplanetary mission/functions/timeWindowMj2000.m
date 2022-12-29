function [dates_MJD2000] = timeWindowMj2000(Trajectory, delta_t)
%
% PROTOTYPE:
%  [dates_MJD2000] = timeWindowMj2000(Trajectory, delta_t);
%
% DESCRIPTION:
%   This function generates a vector of MJD2000 dates between the
%   two defined by the input Trajectory.
%
%  INPUT :
%	Trajectory[2x6] Matrix which contains two dates in the form of a row 
%                   vector defining the earliest and the latest possible
%                   dates.           [year month day hours minutes seconds]
%   delta_t[1]      Variable of type duration, defining the number of days
%                   between two dates.
%
%  OUTPUT:
%	dates_MJD2000[1x?]  Vector containing all the dates between the ones
%	described by the input Trajectory, with a spacing of delta_t days.     [MJD2000]
%
%  FUNCTIONS CALLED:
%   mjd20002date.m
%
% AUTHOR:
%   Virginia Di Biagio Missaglia, Roberto Pistone Nascone, Giuseppe
%   Brentino, Nicolò Galletta 2022
%
% -------------------------------------------------------------------------

%% Convert the input dates in datetime format
t_e = datetime(Trajectory.earliest);
t_l = datetime(Trajectory.latest);

%% Create the vector of dates in datetime forma
dates = t_e : delta_t : t_l;

%% Convert the dates in their original format
dates = datevec(dates);
dates = [dates; Trajectory.latest] ;

%% Convert the dates in MJD2000 format
dates_MJD2000 = zeros (size(dates, 1), 1);
for i = 1 : length(dates_MJD2000)
    dates_MJD2000(i) = date2mjd2000(dates(i, :));
end









