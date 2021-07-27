function T = intervalTo1SecBinTimetable(Start,End,p) 

% intervalToBinTimetable.m

% Simple script to get bin times of specified time from given intervals

Startnum = datenum(Start);
timevec = datevec(Startnum);
tbin = datetime(timevec);
binStartEffort = dateshift(tbin, 'start','second');

Endnum = datenum(End);
timevec = datevec(Endnum);
tbin = datetime(timevec);
binEndEffort = dateshift(tbin, 'start','second');

tbin = [];
effortBin = [];
effortSec = [];
for i = 1: length(Start)
    currInterval = (binStartEffort(i):seconds(1):binEndEffort(i))';
    if ~isempty(currInterval)
    sec = ones(length(currInterval),1);
    % check if first and last bin are complete
    if (binStartEffort(i)-Start(i)) ~= 0
       sec(1) = seconds(seconds(1)-(Start(i)-binStartEffort(i)));
    end
    
    if End(i) > binEndEffort(i)
       sec(end) = seconds(End(i) - binEndEffort(i));
    end
    tbin = [tbin;currInterval];
    effortSec = [effortSec;sec];
    
    end
end

T = timetable(tbin,effortSec);