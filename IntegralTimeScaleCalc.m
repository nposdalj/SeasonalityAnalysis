function its = IntegralTimeScaleCalc(ts)

%Calculate intergral time scale for a time series ts consisting of daily
%whale call counts (in reality, it can be any kind of continuous time
%series, but must be uninterrupted!)
%Note that final output will be in unit increments of your original ts
%bin-size

%based on ITS description in "Data Analysis Methods in PO by Emery & Thomson
%AS 2004

%calculate intergral time scale
avt = mean(ts);
dcc = ts-avt;

[y,lag] = xcorr(dcc,'coeff');
figure(1)
plot(lag,y)

i = ceil(length(y)/2);
its = 0;
while y(i)>0
   its = its+y(i);
   i = i+1;
end

sprintf('ITS is %0.5g days',its)
%note that unit ("days") is based on whatever is your original bin size
disp('Done!')
end