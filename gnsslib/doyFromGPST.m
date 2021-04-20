function doy = doyFromGPST(gpst)
%Approximated day-of-year from GPS-week and tow
gpsWeek=gpst(1);
tow=gpst(2); %(seconds)
JD=gpsWeek*7+2444244.5;
year=(JD-1720981.5)/365.25;
doy=(year-floor(year))*365.25+tow/(3600*24)-30.6001;
if doy<60
    doy=doy-31;
end
end

