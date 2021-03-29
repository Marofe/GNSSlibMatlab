function [stationPos,obs,interval]=readRinexObs(filePath)
%% This function opens RINEX 3.02 observation files.
fileBuffer = fileread(filePath);
ih=strfind(fileBuffer,'END OF HEADER')+length('END OF HEADER');
header=fileBuffer(1:ih);
data=fileBuffer(ih+1:end);
clear fileBuffer
%Initialzie values
interval = -1;
%% Read header
k=0;
lines=strsplit(header,'\n');
while (true)
    k=k+1;
    line = lines{k};                                                   %get line
    splitLine = strsplit(line);                                             %Line splited by spaces
    
    if contains(line,'APPROX POSITION XYZ')                                  % Receiver aprox position
        stationPos=real(str2doubleq(strsplit(line(1:60))));
        stationPos=stationPos(2:end-1)';
        
    elseif contains(line,'SYS / # / OBS TYPES')                    % Observation types for the different constellations (C1C, D1 and S1 only  )
        constellation = line(1);
        if constellation        == 'G'
            hasGps = 1;
            nObservables = real(str2doubleq(line(2:7)));                                  % Number of observables
            observablesGps = splitLine(3:3+nObservables-1);
        elseif constellation    == 'R'
            hasGlonass = 1;
            nObservables = real(str2doubleq(line(2:7)));                                  % Number of observables
            observablesGlonass = splitLine(3:3+nObservables-1);
        elseif constellation    == 'C'
            hasBeidou = 1;
        end
        
    elseif contains(line,'INTERVAL')
        interval=real(str2doubleq(line(5:10)));                       % Measurement intervals (Default 1)
    elseif contains(line,'TIME OF FIRST OBS')
        year = real(str2doubleq(splitLine(2)));
        month = real(str2doubleq(splitLine(3)));
        day = real(str2doubleq(splitLine(4)));
        hour = real(str2doubleq(splitLine(5)));
        minute = real(str2doubleq(splitLine(6)));
        second = real(str2doubleq(splitLine(7)));
        time = [year, month, day, hour, minute, second];
        ti=Date2GPSTime(real(time)); %Transform date to seconds of week
    elseif contains(line,'TIME OF LAST OBS')
        year = real(str2doubleq(splitLine(2)));
        month = real(str2doubleq(splitLine(3)));
        day = real(str2doubleq(splitLine(4)));
        hour = real(str2doubleq(splitLine(5)));
        minute = real(str2doubleq(splitLine(6)));
        second = real(str2doubleq(splitLine(7)));
        time = [year, month, day, hour, minute, second];
        tf=Date2GPSTime(real(time)); %Transform date to seconds of week
    elseif contains(line,'END OF HEADER')
        break;                                                              % End of header loop
    end
end

%default interval=1s
if interval == -1                                               %If itnerval not set interval = 1
    interval = 1;
end
%% Read body
buffer=char(strsplit(data,'\n'));
idEpoch=buffer(:,1)=='>';
nEpochs=nnz(idEpoch);
idEpoch=find(idEpoch);
epochHeader=zeros(nEpochs,4);
year=real(str2doubleq(cellstr(buffer(idEpoch,3:6))));
month=real(str2doubleq(cellstr(buffer(idEpoch,8:9))));
day=real(str2doubleq(cellstr(buffer(idEpoch,11:12))));
hour=real(str2doubleq(cellstr(buffer(idEpoch,14:15))));
min=real(str2doubleq(cellstr(buffer(idEpoch,17:18))));
second=real(str2doubleq(cellstr(buffer(idEpoch,20:29))));
[tow,gpsWeek]=Date2GPSTime(year,month,day,hour,min,second);
epochHeader(:,1)=gpsWeek;
epochHeader(:,2)=tow;
epochHeader(:,3)=real(str2doubleq(cellstr(buffer(idEpoch,32))));
epochHeader(:,4)=real(str2doubleq(cellstr(buffer(idEpoch,33:35))));
epoch=0;
nObs=sum(epochHeader(:,4));
obs=zeros(nObs,12);
fprintf("Reading observation messages...\n");
m=0;
while(epoch<nEpochs-1)
    epoch=epoch+1;
    if ~mod(epoch,floor(nEpochs*0.1))
        fprintf('*');
    end
    data=buffer(idEpoch(epoch)+1:idEpoch(epoch+1)-1,:);
    nSat=epochHeader(epoch,end);
    %     for i = 1:real(nSatellites)                                       % Read the epoch satellites
    %         line = buffer(idEpoch(epoch)+i,:);
    %         if line(1)=='G'
    %             satConst=0;
    %         else
    %             satConst=1;
    %         end
    satConst = data(:,1);
    satConst = satConst=='G'; % Satellites PRN number
    satID = real(str2doubleq(cellstr(data(:,2:3))));                              % Satellites PRN number
    pseudorange = real(str2doubleq(cellstr(data(:,6:17))));
    rangeLLI =real(str2doubleq(cellstr(data(:,18))));
    rangeStrength=real(str2doubleq(cellstr(data(:,19))));
    phase=real(str2doubleq(cellstr(data(:,21:33))));
    phaseLLI = real(str2doubleq(cellstr(data(:,34))));
    phaseStrength=real(str2doubleq(cellstr(data(:,35))));
    doppler=real(str2doubleq(cellstr(data(:,41:49))));
    SNR=real(str2doubleq(cellstr(data(:,60:end))));
    obs(m+1:m+nSat,:) = [gpsWeek(epoch)*ones(nSat,1),tow(epoch)*ones(nSat,1),satConst,satID,pseudorange,rangeLLI,rangeStrength,phase,phaseLLI,phaseStrength,doppler,SNR];   % store data
    m=m+nSat;
end
last=find(obs(:,1)==0,1);
obs=obs(1:last-1,:);
fprintf('\nObservables loaded!\n');