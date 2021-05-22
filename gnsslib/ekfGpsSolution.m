function [time,hx,ns,dop,atmDelay]=ekfGpsSolution(nav,obs,p0,elevMask,atmParam)
% Generate the standard precision solution (~5m) for GPS constellation
% using Extended Kalman Filter
% input:
%       nav         -> broadcast nav messages
%       obs         -> observables (gps only)
%       p0          -> approximated initial position (ECEF)
%       elevMask    -> minimal elevation angle
%% Satellite Orbit
% keplerArray(1)  = svprn;
% keplerArray(2)  = af2;
% keplerArray(3)  = M0;
% keplerArray(4)  = roota;
% keplerArray(5)  = deltan;
% keplerArray(6)  = ecc;
% keplerArray(7)  = omega;
% keplerArray(8)  = cuc;
% keplerArray(9)  = cus;
% keplerArray(10) = crc;
% keplerArray(11) = crs;
% keplerArray(12) = i0;
% keplerArray(13) = idot;
% keplerArray(14) = cic;
% keplerArray(15) = cis;
% keplerArray(16) = Omega0;
% keplerArray(17) = Omegadot;
% keplerArray(18) = toe;
% keplerArray(19) = af0;
% keplerArray(20) = af1;
% keplerArray(21) = week_toe;
% keplerArray(22) = tgd;
% keplerArray(23) = txTime;
% keplerArray(24) = toc;
% obs=[gpsWeek,tow,satConst,satID,pseudorange,rangeLLI,rangeStrength,phase,phaseLLI,phaseStrength,doppler,SNR];
%%
fprintf('Processing gnss observations...\n')
time=unique(obs(:,2));
Ns=length(time); %Number of GNSS samples
dop=zeros(Ns,4);
ns=zeros(Ns,1);
atmDelay=zeros(Ns,2);
hx=zeros(Ns,4);
P=zeros(4,4,Ns);
hx(1,:)=[p0' 0];
P(:,:,1)=diag([10 10 10 1e6]);
A=eye(4);
Q=eye(4)*100;
for k=1:Ns-1
    %% Epoch
    epoch=obs(obs(:,2)==time(k),:);
    %% Update
    [hx(k,:),P(:,:,k),ns(k),ionT,tropT,lla]=iteratedUpdate(hx(k,:)',P(:,:,k),epoch,nav,atmParam,time(k),elevMask);
    %% Prediction
    hx(k+1,:)=hx(k,:)*A';
    P(:,:,k+1)=A*P(:,:,k)*A'+Q;
    %% Compute the Dilusion of Precision (DOP)
    Cen=DCM_en(lla(1),lla(2));
    Pn=Cen'*P(1:3,1:3,k)*Cen; % ECEF->NED
    dop(k,1)=sqrt(trace(Pn)); %Geometric (GDOP)
    dop(k,2)=sqrt(Pn(1,1)+Pn(2,2)); %Horizontal (HDOP)
    dop(k,3)=sqrt(Pn(3,3)); %Vertical (VDOP)
    dop(k,4)=sqrt(P(4,4,k)); %Time (TDOP)
    %% Atmospheric delay
    atmDelay(k,1)=tropT/ns(k); %ZTD (Zenith Tropospheric Delay)
    atmDelay(k,2)=ionT/ns(k); %Ionospheic
    %
    if ~mod(k,round(Ns/10))
        fprintf('*')
    end
end
fprintf('\nProcessing finished!')
end