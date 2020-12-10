function T = intervalToBinTimetable_Only1RowEffort(Start,End,p) 


Startnum = datenum(Start);
timevec = datevec(Startnum);
[m,n] = size(timevec);
binStartEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(m)]);

Endnum = datenum(End);
timevec = datevec(Endnum);
[m,n] = size(timevec);
binEndEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(m)]);

tbin = [];
for i = 1: length(Start)
    tbin = [tbin;(binStartEffort(i):minutes(5):binEndEffort(i))']; 
end

T = timetable(tbin);
T.bin = ones(height(T),1);