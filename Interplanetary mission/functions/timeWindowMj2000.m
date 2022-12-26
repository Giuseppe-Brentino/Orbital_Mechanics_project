function [dates_MJD2000] = timeWindowMj2000(Trajectory, delta_t)
%{
 This function generates the vector of MJD2000 dates between the earliest
 and latest defined in the struct Trajectory.
 
 INPUTS:
 Trajectory [2 x 6]: Matrix which contains two dates in the form of a row 
    vector [year month day hours minutes seconds] defining the earliest and the
    latest possible dates
 delta_t [1]: variable of type duration defining the number of days between 
    two dates.

 OUTPUTS:
 dates_MJD2000 [1x?]: vector containing al the dates (in MJD2000 form) between the
    ones described by the input Trajectory, with a spacing of delta_t days
 
 Contributors:
 Giuseppe Brentino

%}

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









