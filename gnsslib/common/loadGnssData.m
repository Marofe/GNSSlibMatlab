function [gnss,ti,tf,Fs] = loadGnssData(fileName)
%% Load data from RTKLIB Post-Processing
data = load(fileName);
    %% RTKLIB FILE .pos
    % ECEF position
    week=data(:,1);
    time=data(:,2); %seconds
    x=data(:,3); %meters
    y=data(:,4); %meters
    z=data(:,5); %meters
    fix=data(:,6); %1=fix; 2=float
    ns=data(:,7);
    sdv=data(:,8:13);
    age=data(:,14);
    ratio=data(:,15);
    % GPS time (GPST)
    dz=abs(diff(z));
    Ni=find(dz>1,1);
    Nf=length(dz)-find(flip(dz)>1,1);
    ti=time(Ni)-300;
    tf=time(Nf)+300;  
    Fs=1/round(mean(diff(time)));
    % Transform ECEF to Geodetic
    [lat, lon, alt]=llaFromEcef(x,y,z);
    Gaps=sum(diff(time)>1.5*median(diff(time)))-1;
    if Gaps<0
        Gaps=0;
    end
    fprintf("GNSS GAPs: %d (%.2f%%)\n",Gaps,Gaps/size(data,1)*100);
    fprintf('GNSS Quality: %.2f%%\n',sum(fix==1)/length(time)*100)
    gnss=[time x y z lat lon alt fix ns sdv age ratio];
end

