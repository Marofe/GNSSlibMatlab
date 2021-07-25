function [orbit]=readSP3file(filePath)
%% This function opens SP3-c precise orbit files.
fileBuffer = fileread(filePath);
lines=strsplit(fileBuffer,'\n');
clear fileBuffer
header=char(lines{1:22});
data=char(lines{23:end});
%% Read header

%% Read body
idEpoch=data(:,1)=='*';
nEpochs=nnz(idEpoch);
idEpoch=find(idEpoch);
epochHeader=zeros(nEpochs,2);
year=real(str2doubleq(cellstr(data(idEpoch,4:7))));
month=real(str2doubleq(cellstr(data(idEpoch,9:10))));
day=real(str2doubleq(cellstr(data(idEpoch,12:13))));
hour=real(str2doubleq(cellstr(data(idEpoch,15:16))));
min=real(str2doubleq(cellstr(data(idEpoch,18:19))));
second=real(str2doubleq(cellstr(data(idEpoch,21:31))));
[tow,gpsWeek]=Date2GPSTime(year,month,day,hour,min,second);
epochHeader(:,1)=gpsWeek;
epochHeader(:,2)=tow;
epoch=0;
obs=zeros(nEpochs,8);
fprintf("Reading SP3 file...\n");
for i=1:32
buffer=data(idEpoch+i,:);
x=real(str2doubleq(cellstr(buffer(:,5:18)))); %km
y=real(str2doubleq(cellstr(buffer(:,19:32)))); %km
z=real(str2doubleq(cellstr(buffer(:,33:46)))); %km
t=real(str2doubleq(cellstr(buffer(:,47:60)))); %micro-s
sdx=real(str2doubleq(cellstr(buffer(:,62:63)))); %mm
sdy=real(str2doubleq(cellstr(buffer(:,65:66)))); %mm
sdz=real(str2doubleq(cellstr(buffer(:,68:69)))); %mm
sdt=real(str2doubleq(cellstr(buffer(:,71:73)))); %pico-s
orbit(:,:,i)=[tow [x y z]*1e3 t sdx sdy sdz sdt];
end
fprintf('\nPrecise orbit file loaded!\n');