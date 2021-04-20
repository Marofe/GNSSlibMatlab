function [gnss,ti,tf,Fs] = loadGlabGnssData(fileName)
%% Load data from gLab Post-Processing
fileBuffer = fileread(fileName);
% lines=char(strsplit(fileBuffer,'\n'));
lines=strsplit(fileBuffer,'\n');
N=numel(lines);
x=zeros(N,1);
y=zeros(N,1);
z=zeros(N,1);
for k=1:N
    line=strsplit(lines{k});
    x(k)=real(str2doubleq(line{6}));%data(:,6); %meters
    y(k)=real(str2doubleq(line{7}));%data(:,7); %meters
    z(k)=real(str2doubleq(line{8}));%data(:,8); %meters
    %     lat=data(:,15);
    %     lon=data(:,16);
    %     alt=data(:,17);
    %     ns=data(:,32);
    %     sdv=data(:,12:14);
    %     % GPS time (GPST)
    %     time=1:size(data,1);
end
[lat, lon, alt]=llaFromEcef(x,y,z);
gnss=[x y z lat lon alt];
end

