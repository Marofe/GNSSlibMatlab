function [time,x,ns,dop,sdNa,baseline]=KinematicFloatSolution(ephemeris,obsBase,obsRover,pbase,p0rover,elevMask)
% Generate the high precision solution (~5cm) for GPS constellation
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
%% Non-zeros phase/range measurement
obsRover(obsRover(:,8)==0,:)=[];
obsBase(obsBase(:,8)==0,:)=[];
obsRover(obsRover(:,5)==0,:)=[];
obsBase(obsBase(:,5)==0,:)=[];
%% Nav data
nav=ephemeris.gpsEphemeris;
atmParam=ephemeris.ionosphericParameters;
%%
fprintf('\nProcessing Kinematic Mode...\n')
time=unique(obsRover(:,2));
N=length(time);
sats=unique(nav(1,:));
Nsats=numel(sats);
%% Allocate memmory
x=zeros(6,N);
xb=zeros(6,N);
Px=zeros(6,6,N);
dop=zeros(N,3);
ns=zeros(N,1);
baseline=zeros(N,1);
sdNa=zeros(Nsats,N); %single-difference ambiguity
satsOnView=zeros(Nsats,N); %sats on
sdNab=zeros(Nsats,N); %single-difference ambiguity
Pna=zeros(Nsats,N); %single-difference ambiguity variance
%%
iterMax=15;
%% Dynamic model (Const. Velocity)
dt=1; %1hz GPS
A=[eye(3) dt*eye(3);zeros(3) eye(3)];
Q=blkdiag(0*eye(3),.1*eye(3)); 
%% initial position/velocity
x(:,1)=[p0rover;zeros(3,1)];
Px(:,:,1)=10*eye(6);
Pna(:,1)=100*ones(Nsats,1);
for k=1:N
    %% Epochs
    trcv=time(k); %time of epoch (toe)
    epoch=obsRover(obsRover(:,2)==trcv,:);
    epochBase=obsBase(obsBase(:,2)==trcv,:);
    iter=1;
    err=inf;
    hatp=x(1:3,k); %prior
    hatv=x(4:6,k);
    Na0=sdNa(:,k);
    xx0=[hatp;hatv]; %prior
    while err>1e-3 && iter<=iterMax
        lla=SingleLlaFromEcef(hatp);
        z1=[];
        z2=[];
        satsOn=[];
        acceptSatRef=0;
        while ~acceptSatRef
            satId1=epoch(1,4); %satref (first of the epoch)
            base1=epochBase(epochBase(:,4)==satId1,:);
            if size(base1,1)==0
                %satref not found in the base obs
                epoch=epoch(2:end,:);
            else
                %satref found in the base obs
                acceptSatRef=1;
            end
        end
        
        %rover ref measurements
        rho1_r=epoch(1,5); %pseudo-range sat1 (rover)
        Phi1_r=epoch(1,8)*lambda1; %phase-range sat1 (rover)
        SNR1_r=epoch(1,end); %SNR sat1 (rover)
        
        %compute satref position
        id=find(nav(1,:)==satId1);  %select all ephemeris for the ref satellite
        [~,j]=min(abs(trcv-nav(18,id)));  %seek for the most recent orbit parameters
        j=id(j);
        tems=trcv-rho1_r/c; %emission time (GPST)
        psat1=satPosition(nav(:,j),tems);
        theta=wie*norm(psat1-hatp)/c;
        psat1=rotZ(theta)*psat1;
        Cen=DCM_en(lla(1),lla(2));
        u1=(psat1-hatp)/norm(psat1-hatp);
        los1=Cen'*u1; %Line of Sight (LOS) (rover)
        %azim=atan2(los(1),los(2)); %Azimute (NED frame)
        elev1=asin(-los1(3)); %Elevation (rad) (NED frame)
        u=u1';
        sigmaR1=10/SNR1_r+3*exp(-2*elev1/(pi/2));
        R=sigmaR1*sigmaR1';
        
        %base ref measurements
        rho1_b=base1(5); %pseudo-range sat1 (base)
        Phi1_b=base1(8)*lambda1; %phase-range sat1 (base)
        
        Ni=size(epoch,1); %number of satellites (=measurements)
        for m=2:Ni %start at 2 because the first is the ref
            satId=epoch(m,4);
            base=epochBase(epochBase(:,4)==satId,:);
            if size(base,1)>0
                %satId found in the base epoch
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
                %% Compute Line-of-Sight (LOS)
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
                    bar_ddrho=(norm(psat1-hatp)-norm(psat1-pbase))-(norm(psat-hatp)-norm(psat-pbase));
                    hat_ddrho=bar_ddrho; %predicted pseudo-range
                    hat_ddPhi=bar_ddrho+lambda1*(sdNa(sats==satId1,k)-sdNa(sats==satId,k));  %predicted phase-range
                    %% Residuals
                    z1=[z1;ddPhi-hat_ddPhi]; %phase-range residual
                    z2=[z2;ddrho-hat_ddrho]; %pseudo-range residual
                    satsOn=[satsOn satId];
                    u=[u;range'/bar_rho];
                    sigmaR=10/SNR_r+3*exp(-2*elev/(pi/2));
                    R=blkdiag(R,sigmaR*sigmaR');
                end
            end
        end
        %% Weighted-Least-Squared-Error Estimation
        pref=hatp;
        Nii=numel(satsOn); %(2...M)
        D=[ones(Nii,1) -eye(Nii)];
        H=[-D*u zeros(Nii,3) lambda1*D;-D*u zeros(Nii,Nii+1+3)];
        z=[z1;z2];
        R=blkdiag(D*R/1000*D',D*R/10*D');
        hatx=[hatp;hatv;zeros(Nii+1,1)]; %ith-prior
        x0=[xx0;zeros(Nii+1,1)]; %first-prior
        satsOn=[satId1 satsOn];
%        z3=zeros(Nii+1,1);
        Pnai=zeros(Nii+1,1);
        
        for s=1:Nii+1
                hatx(6+s)=sdNa(sats==satsOn(s),k);
                x0(6+s)=Na0(sats==satsOn(s));
                Pnai(s)=Pna(sats==satsOn(s),k);
                satsOnView(sats==satsOn(s),k)=1;
        end
        
        P0=blkdiag(Px(:,:,k),diag(Pnai));
        K=P0*H'/(R+H*P0*H');
        hatx=x0+K*(z+H*(hatx-x0));
        hatp=hatx(1:3);
        hatv=hatx(4:6);
        err=norm(hatp-pref);
        if err>1e3
            %outlier
            hatp=pref;
        else
            for s=1:Nii+1
              sdNa(sats==satsOn(s),k)=hatx(6+s);
            end
        end
        iter=iter+1;
    end
    %% result
    x(:,k)=hatx(1:6);
    baseline(k)=norm(pbase-x(1:3,k));
    P=P0-K*(R+H*P0*H')*K';
    Px(:,:,k)=P(1:6,1:6);
    Pnai=P(7:end,7:end);
    ns(k)=size(u,1);
    %% Ambiguity resolution
%     hatN=D*hatx(7:end);
%     G=blkdiag(eye(6),D);
%     W=G*P*G';
%     WN=W(7:end,7:end);
%     WXN=W(1:6,7:end);
%     optN=round(hatN);
%     x(:,k)=x(:,k)+WXN/WN*(optN-hatN);
%     for s=1:Nii+1
%         sdNa(sats==satsOn(s),k)=optN(s);
%     end
    %% Compute the Dilusion of Precision (DOP)
    Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
    dop(k,1)=sqrt(trace(Pn)); %Geometric (GDOP)
    dop(k,2)=sqrt(Pn(1,1)+Pn(2,2)); %Horizontal (HDOP)
    dop(k,3)=sqrt(Pn(3,3)); %Vertical (VDOP)
    if ~mod(k,round(length(time)/10))
        fprintf('*')
    end
    %% time-update (prediction)
    if k<N
        x(:,k+1)=A*x(:,k);
        sdNa(:,k+1)=sdNa(:,k);
        Px(:,:,k+1)=A*Px(:,:,k)*A'+Q;
        Pnai=Pnai+.1*eye(Nii+1);
        Pnai=diag(Pnai);
        for s=1:Nii+1
            Pna(sats==satsOn(s),k+1)=Pnai(s);
        end
    end
end
%% Smoothing
% xb(:,end)=x(:,end);
% for k=N-1:-1:1
%     Nii=sum(satsOnView(:,k)==1);
%     hx=[x(:,k);zeros(Nii,1)]; %hx(k)
%     for s=1:Nii+1
%         Pna(sats==satsOn(s),k+1)=Pnai(s);
%     end
% hx=A
% end
%% solution
x=x';
fprintf('\nKinematic processing finished!\n')
fprintf('Max. baseline: %.2fkm\n',max(baseline)/1e3)
if max(baseline)/1e3 >10
    warning('Long-baseline!')
end
end