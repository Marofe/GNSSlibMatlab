close all
clear all
clc
format long
%%
fileId = fopen('neo6m_gnss_log.txt','r');
%$GPGGA,195646.00,2249.62941,S,04704.69851,W,1,10,1.02,606.0,M,-5.7,M,,*46
i=0;
line=fgetl(fileId);
while ischar(line)
line=fgetl(fileId);
if length(line)>2 && contains(line,'$GPGGA')
    raw=split(line,',');
    i=i+1;
    gnss(i,1)=datenum(raw(2),'HHMMSS'); %GPST
    gnss(i,2)=-str2double(raw(3))/100; %LAT
    gnss(i,3)=-str2double(raw(5))/100; %LON
    gnss(i,4)=str2double(raw(10)); %ALT
    deg=floor(abs(gnss(i,2:3))).*sign(gnss(i,2:3));
    gnss(i,2:3)=deg+(gnss(i,2:3)-deg)/60*100;
end
end
fclose(fileId);
[x,y,z]=ecefFromLLA(gnss(:,2),gnss(:,3),gnss(:,4));
[N,E,D]=nedFromEcef(x,y,z,gnss(1,2),gnss(1,3));
gnss=[gnss N' E' D'];
%%
figure
plot(gnss(:,3),gnss(:,2),'ob')
grid on
xlabel('lon (deg)')
ylabel('lat (deg)')
title('Geographic Position')
figure
hold on
plot(gnss(:,5)-mean(gnss(~isnan(gnss(:,5)),5)),gnss(:,6)-mean(gnss(~isnan(gnss(:,6)),6)),'ro')
stdX=std(gnss(~isnan(gnss(:,5)),2))
stdY=std(gnss(~isnan(gnss(:,6)),3))
xlabel('north (m)')
ylabel('east (m)')
grid on
title('NED Position')
figure
plot(gnss(:,1),gnss(:,4),'ro')
zm=mean(gnss(~isnan(gnss(:,4)),4));
hold on
plot([gnss(1,1) gnss(end,1)],[zm zm],'g','linewidth',2)
title('Altitude')
xlabel('GPST')
ylabel('alt (m)')
grid on