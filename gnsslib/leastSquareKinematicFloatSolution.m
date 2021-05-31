function [time,p,ns,dop,sdNa]=leastSquareKinematicFloatSolution(nav,obsBase,obsRover,pbase,p0rover,elevMask,atmParam)
% Generate the high precision solution (~5m) for GPS constellation
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
fprintf('\nProcessing Kinematic Mode...\n')
time=unique(obsRover(:,2));
N=length(time);
sats=unique(nav(1,:));
Nsats=numel(sats);
%% Allocate memmory
p=zeros(N,3);
dop=zeros(N,3);
ns=zeros(N,1);
sdNa=zeros(Nsats,N+1);
%%
iterMax=15;
for k=1:N
    %% Epoch
    trcv=time(k);
    epoch=obsRover(obsRover(:,2)==trcv,:);
    epochBase=obsBase(obsBase(:,2)==trcv,:);
    if k>1
        hatp=p(k-1,:)'; %hatp=[hxr,hyr,hzr]
    else
        hatp=p0rover;
    end
    iter=1;
    err=inf;
    while err>1e-3 && iter<=iterMax
        lla=SingleLlaFromEcef(hatp);
        z1=[];
        z2=[];
        G=[];
        W=[];
        %ref sat = the one with highest elevation
        satsOn=[];
        satId1=epoch(1,4);
        base1=epochBase(epochBase(:,4)==satId1,:);
        if size(base1,1)==0
            epoch=epoch(2:end,:);
        end
        satId1=epoch(1,4);
        base1=epochBase(epochBase(:,4)==satId1,:);
        rho1_r=epoch(1,5); %pseudo-range sat1 (rover)
        Phi1_r=epoch(1,8)*lambda1; %phase-range sat1 (rover)
        id=find(nav(1,:)==satId1);  %select all ephemeris for the satellite
        [~,j]=min(abs(trcv-nav(18,id)));  %seek for the most recent orbit parameters
        j=id(j);
        tems=trcv-rho1_r/c; %(GPST)
        psat1=satPosition(nav(:,j),tems);
        theta=wie*norm(psat1-pbase)/c;
        psat1=rotZ(theta)*psat1;
        rho1_b=base1(5); %pseudo-range sat1 (base)
        Phi1_b=base1(8)*lambda1; %phase-range sat1 (base)
        Ni=size(epoch,1); %number of satellites (=measurements)
        for m=2:Ni
            satId=epoch(m,4);
            base=epochBase(epochBase(:,4)==satId,:);
            if size(base,1)>0
            id=find(nav(1,:)==satId);  %select all ephemeris for the satellite
            [~,j]=min(abs(time(k)-nav(18,id)));  %seek for the most recent orbit parameters
            j=id(j);
            %toe=nav(18,j); %time of ephemeris
            %% Observables
            %rover
            rho_r=epoch(m,5); %pseudo-range (~5m) (rover)
            Phi_r=epoch(m,8)*lambda1; %phase-range (~10cm+ambiguity) (rover)
            %doppler=-epoch(m,11)*lambda1; %range-rate
            SNR_r=epoch(m,end); %Signal-to-Noise ratio (rover)
            
            %base
            rho_b=base(5); %pseudo-range (~5m) (base)
            Phi_b=base(8)*lambda1; %phase-range (~10cm+ambiguity) (base)
            %doppler=-epoch(m,11)*lambda1; %range-rate
            SNR_b=base(end); %Signal-to-Noise ratio (base)
            
            %% Compute the emission time (approx)
            tems=trcv-rho_r/c; %(GPST)
            
            %% Compute satellite position (at emission time)
            % From Broadcast nav msg
            [psat0,E,a]=satPosition(nav(:,j),tems);
            range0=psat0-hatp;
            bar_rho0=norm(range0); %approx. geometric range
            
            %% Take in account the Earth's rotation
            theta=wie*bar_rho0/c;
            psat=rotZ(theta)*psat0;
            range=psat-hatp;
            bar_rho=norm(range); %approx. geometric range
            Cen=DCM_en(lla(1),lla(2));
            los=Cen'*range/bar_rho; %Line of Sight (LOS)
            %azim=atan2(los(1),los(2)); %Azimute (NED frame)
            elev=asin(-los(3)); %Elevation (rad) (NED frame)
            elevd=rad2deg(elev); %Elevation (deg)
             if elevd>elevMask
                %% Compute Double Difference
                ddrho=(rho1_r-rho1_b)-(rho_r-rho_b); %rho_rb_1j
                ddPhi=(Phi1_r-Phi1_b)-(Phi_r-Phi_b); %rho_rb_1j
                %% Pseudo-range
                bar_ddrho=norm(psat1-hatp)-norm(psat-hatp)-(norm(psat1-pbase)-norm(psat-pbase));
                hat_ddrho=bar_ddrho; %predicted pseudo-range
                hat_ddPhi=bar_ddrho+lambda1*sdNa(sats==satId,k);
                %% Residuals
                z1=[z1;ddPhi-hat_ddPhi]; %phase-range residual
                z2=[z2;ddrho-hat_ddrho]; %pseudo-range residual
                satsOn=[satsOn satId];
                r1=-(psat1-hatp)/norm(psat1-hatp);
                r=-(psat-hatp)/norm(psat-hatp);
                G=[G;r1'-r'];
                sigmaR=10/SNR_r+3*exp(-2*elev/(pi/2));
                R=sigmaR*sigmaR';
%                 a=0.13;b=0.53;c=10; %from MOPS
%                 sigmaR=a+b*exp(-rad2deg(elev)/c);
                W=blkdiag(W,1/R);
             end
            end
        end
        %% Weighted-Least-Squared-Error Estimation
        pref=hatp;
        Nii=numel(satsOn);
        D=[ones(Nii,1) -eye(Nii)];
        H=[G lambda1*D;G zeros(Nii,Nii+1)];
        z=[z1;z2];
        W=blkdiag(W*100,W);
        hatx=[hatp;zeros(Nii+1,1)];
        if iter>1
            for s=1:Nii+1
                hatx(3+s)=sdNa(sats==satsOn(s),k);
            end
        end
        iP=H'*W*H;
        iP=0.5*(iP+iP');
        hatx=hatx+iP\H'*W*z;
        hatp=hatx(1:3);
        for s=1:Nii+1
            sdNa(sats==satsOn(s),k+1)=hatx(3+s);
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
    ns(k)=size(G,1);
    P=inv(H'*W*H);
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
fprintf('\n Kinematic processing finished!')
end