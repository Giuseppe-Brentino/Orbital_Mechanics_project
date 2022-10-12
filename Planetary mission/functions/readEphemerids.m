function [kep] = readEphemerids(filename)

ephemerids=readmatrix(filename);

found=0;
row=3;

while ~found
    if isnan(ephemerids(row,1)) 
        found=1;
    end
    row=row+1;
end

row=row-2;

%[a,e,i,OM,om,theta]
kep=[ephemerids(3:row,12), ephemerids(3:row,3), ...
    ephemerids(3:row,5:7), ephemerids(3:row,11)];
