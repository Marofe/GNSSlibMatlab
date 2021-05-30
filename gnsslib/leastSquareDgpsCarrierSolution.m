function [time,p,ns,dop]=leastSquareDgpsCarrierSolution(nav,obsBase,obsRover,p0base,p0rover,elevMask,atmParam)
% Generate the standard precision solution (~5m) for GPS constellation
% input:
%       nav         -> broadcast nav messages
%       obs         -> observables (gps only)
%       p0          -> approximated initial position (ECEF)
%       elevMask    -> minimal elevation angle
%% Constants
% mu=3986004.418e8;% gravitational constant (m^3/s^2)
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
f1=1.57542*1e9; %L1 GPS frequency
lambda1=c/f1; %L1 GPS wavelength
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
fprintf('\nProcessing DGPS-Carrier...\n')
time=unique(obsRover(:,2));
N=length(time);
sats=unique(nav(1,:));
Nsats=numel(sats);
%% Allocate memmory
p=zeros(N,4);
dop=zeros(N,4);
ns=zeros(N,1);
Na=zeros(Nsats,1);
%%
iterMax=15;
for k=1:N
    %% Epoch
    trcv=time(k);
    epoch=obsRover(obsRover(:,2)==trcv,:);
    epochBase=obsBase(obsBase(:,2)==trcv,:);
    Ni=size(epoch,1);
    if k>1
        cdt=p(k-1,4); %relative clock-offset (rover-base)
        hatp=p(k-1,:)';
    else
        cdt=0; %relative clock-offset (rover-base)
        hatp=[p0rover;cdt];
    end
    iter=1;
    pref=zeros(4,1);
    err=inf;
    while err>1e-3 && iter<=iterMax
        lla=SingleLlaFromEcef(hatp(1:3));
        res1=[];
        res2=[];
        H=[];
        W=[];
        satsOn=[];
        for m=1:Ni
            satId=epoch(m,4);
            base=epochBase(epochBase(:,4)==satId,:);
            id=find(nav(1,:)==satId);  %select all ephemeris for the satellite
            [~,j]=min(abs(time(k)-nav(18,id)));  %seek for the most recent orbit parameters
            j=id(j);
            %toe=nav(18,j); %time of ephemeris
            %TGD=nav(22,j); %Total Group Delay
            %% Observables
            r=epoch(m,5); %pseudo-range (~5m)
            phase=epoch(m,8)*lambda1; %phase-range (~10cm+ambiguity)
            %doppler=-epoch(m,11)*lambda1; %range-rate
            SNR=epoch(m,end); %Signal-to-Noise ratio
            %% Compute the emission time (approx)
            tems=trcv-r/c; %(GPST)
            %% Compute satellite position (at emission time)
            % From Broadcast nav msg
            [satPos0,E,a]=satPosition(nav(:,j),tems);
            range0=satPos0-hatp(1:3);
            rho0=norm(range0); %approx. geometric range
            %% Take in account the Earth's rotation
            theta=wie*rho0/c;
            satPos=rotZ(theta)*satPos0;
            %rhoSat=norm(satPos);
            range=satPos-hatp(1:3);
            rho=norm(range); %approx. geometric range
            %rhoRcv=norm(hatp(1:3)); %predicted rcv position
            Cen=DCM_en(lla(1),lla(2));
            los=Cen'*range/rho; %Line of Sight (LOS)
            %azim=atan2(los(1),los(2)); %Azimute (NED frame)
            elev=asin(-los(3)); %Elevation (rad) (NED frame)
            elevd=rad2deg(elev); %Elevation (deg)
             if elevd>elevMask
                %% Compute correction
                r0=norm(p0base-satPos);
                deltaRho=base(5)-r0;
                deltaPhase=base(8)*lambda1-r0;
                %% Pseudo-range
                r=r-deltaRho; %single-differenced pseudo-range
                phase=phase-deltaPhase; %single-differenced phase
                tilR=rho+cdt; %predicted pseudo-range
                tilPhase=rho+cdt+lambda1*Na(sats==satId);
                res1=[res1;r-tilR]; %pseudo-range residual
                res2=[res2;phase-tilPhase]; %phase-range residual
                satsOn=[satsOn satId];
                H=[H;-range'/rho 1];
                sigmaR=10/SNR+3*exp(-2*elev/(pi/2));
                R=sigmaR*sigmaR';
%                 a=0.13;b=0.53;c=10; %from MOPS
%                 sigmaR=a+b*exp(-rad2deg(elev)/c);
                W=blkdiag(W,1/R);
             end
        end
        %% Weighted-Least-Squared-Error Estimation
        pref=hatp;
        Nii=numel(satsOn);
        HH=[H zeros(Nii);H eye(Nii)*lambda1];
        res=[res1;res2];
        W=blkdiag(W,W*100);
        hatx=[hatp;zeros(Nii,1)];
        if iter>1
            for s=1:Nii
                hatx(4+s)=Na(sats==satsOn(s));
            end
        end
        PP=HH'*W*HH;
        PP=0.5*(PP+PP');
        hatx=hatx+PP\HH'*W*res;
        hatp=hatx(1:4);
        cdt=hatp(4);
        for s=1:Nii
            Na(sats==satsOn(s))=hatx(4+s);
        end
        err=norm(hatp-pref);
        iter=iter+1;
    end
    %% Consistency test
%     satValid=abs(res)<2.5;
%     if sum(~satValid)~=0
%         W=diag(W);
%         W=diag(W(satValid));
%         H=H(satValid,:);
%         hatp=hatp+(H'*W*H)\H'*W*res(satValid);
%     end
    %% result
    p(k,:)=hatp';
    ns(k)=size(H,1);
    P=inv(HH'*W*HH);
    %% Compute the Dilusion of Precision (DOP)
    Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
    dop(k,1)=sqrt(trace(Pn)); %Geometric (GDOP)
    dop(k,2)=sqrt(Pn(1,1)+Pn(2,2)); %Horizontal (HDOP)
    dop(k,3)=sqrt(Pn(3,3)); %Vertical (VDOP)
    dop(k,4)=sqrt(P(4,4)); %Time (TDOP)
    if ~mod(k,round(length(time)/10))
        fprintf('*')
    end
end
fprintf('\n DGPS processing finished!')
end