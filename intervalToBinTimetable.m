function T = intervalToBinTimetable(Start,End,p) 


Startnum = datenum(Start);
timevec = datevec(Startnum);
[q,~] = size(timevec);
binStartEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(q,1)]);

Endnum = datenum(End);
timevec = datevec(Endnum);
[q,~] = size(timevec);
binEndEffort = datetime([timevec(:,1:4), floor(timevec(:,5)/p.binDur)*p.binDur, ...
    zeros(q,1)]);

tbin = [];
for i = 1: length(Start)
    tbin = [tbin;(binStartEffort(i):minutes(5):binEndEffort(i))']; 
end

T = timetable(tbin);
T.bin = ones(height(T),1);