function out=KinematicSolution(ephemeris,obsBase,obsRover,pbase,p0rover,param)
% Generate the high precision solution (~5cm) for GPS constellation
% input:
%       nav         -> broadcast nav messages
%       obs         -> observables (gps only)
%       p0          -> approximated initial position (ECEF)
%       elevMask    -> minimal elevation angle
fprintf('Processing Kinematic Mode...\n')
%% constants
c=2.99792458e8; %Speed of light in vacuum (m/s)
f1=1.57542*1e9; %L1 GPS frequency
lambda1=c/f1; %L1 GPS wavelength
%% Clean Non-zeros phase/range measurement && LLI (loss of lock indicator)
obsRover(obsRover(:,8)==0,:)=[];
obsBase(obsBase(:,8)==0,:)=[];
obsRover(obsRover(:,5)==0,:)=[];
obsBase(obsBase(:,5)==0,:)=[];
obsBase(obsBase(:,9)==0,9)=NaN; %phase LLI
obsRover(obsRover(:,6)==0,6)=NaN; %code LLI
%% Nav data
nav=ephemeris.gpsEphemeris;
%%
time=unique(obsRover(:,2)); %GPST (receiver)
N=length(time); %number of epochs
sats=unique(nav(1,:)); %all sats within the dataset
Nsats=numel(sats); %max. number of satellites
%% Allocate memmory
hx=zeros(6,N);
Px=zeros(6,6,N);
hatNa=zeros(Nsats,N); %single-difference ambiguity
Pna=zeros(Nsats,Nsats,N); %single-difference ambiguity covariance
dop=zeros(N,3);
cycleSlip=zeros(N,1);
satsOnView=zeros(Nsats,N); %sats on
satsTime=zeros(Nsats,1); %sats on
refSatId=zeros(N,1); %sats ref (with highest elevation)
Ra=zeros(N,1); %ratio-test
trPna=zeros(N,1); %trace of ambiguity covariance
res_phi=zeros(Nsats,N); %residues phase
res_rho=zeros(Nsats,N); %residues code
%% params
iterMax=5; %maximum iterations of the ILS-update step
dt=1; %1hz GPS
sigmaVn=1; %north-acc noise (m/s/sqrt(s))
sigmaVe=1; %east-acc noise (m/s/sqrt(s))
sigmaVd=5; %down-acc noise (m/s/sqrt(s))
sigmaN=0.0001; %Ambiguity fix-and-hold
sigmaNfloat=.5; %Ambiguity var when float
sigmaNfix=1e-5; %Ambiguity var when fix
Rthres=3; %Ratio-test
%code/phase std = a+b/sin(elev)
a=0.1; %Independent noise
b=0.5; %Elev contribution
%Initial variance
p0=[1000*ones(1,3),100*ones(1,3)]; %position & velocity initial variance
pN0=10000; %initial ambiguity & cycle slip variance
flagFix=1; %1=try to fix; 0=skip
flagN0=1; %reset initial ambiguity to zero
alpha=0.9;
gN=100;
unpackStruct(param);
sigmaNp=sigmaNfloat*ones(Nsats,1); %Initial Ambiguity var
%% Dynamic model (Const. Velocity)
Qvn=diag([sigmaVn^2,sigmaVe^2,sigmaVd^2]);
%% Initialization
hx(:,1)=[p0rover;zeros(3,1)];
hatNa(:,1)=zeros(Nsats,1);
Px(:,:,1)=diag(p0);
Pna(:,:,1)=pN0*eye(Nsats);
%% Iterated-EKF-Processing
for k=1:N
    iter=1;
    err=inf;
    hatp=hx(1:3,k); %prior position
    hatv=hx(4:6,k); %prior velocity
    hatx0=hx(:,k);
    Na0=hatNa(:,k);
    %% Epochs
    %time(k) => time of epoch (toe)
    epoch=obsRover(obsRover(:,2)==time(k),:);
    epochBase=obsBase(obsBase(:,2)==time(k),:);
    % obs=[gpsWeek,tow,satConst,satID,pseudorange,rangeLLI,rangeStrength,phase,phaseLLI,phaseStrength,doppler,SNR];
    %% ILS
    while err>1e-3 && iter<=iterMax
        data=getSatellitesInfo(epoch,epochBase,nav,hatp,pbase,hatNa(:,k),sats,elevMask);
        Ni=size(data,1);
        refSatId(k)=data(1,1);
        lla=SingleLlaFromEcef(hatp);
        Cen=DCM_en(lla(1),lla(2));
        
        dRho=data(:,3);
        dPhi=data(:,4);
        hat_dRho=data(:,13);
        hat_dPhi=data(:,14);
        elev=data(:,2);
        R_dRho=2*diag(a+b./sin(elev));
        R_dPhi=R_dRho/(1e4);
        
        D=[ones(Ni-1,1) -eye(Ni-1)];
        ddRho=D*dRho;
        ddPhi=D*dPhi;
        hat_ddRho=D*hat_dRho;
        hat_ddPhi=D*hat_dPhi;
        z_ddRho=ddRho-hat_ddRho;
        R_ddRho=D*R_dRho*D';
        R_ddPhi=D*R_dPhi*D';
        outlier=abs(z_ddRho)>3*sqrt(diag(R_ddRho));
        z_ddRho(outlier)=[];
        D_rho=D(~outlier,:);
        R_ddRho=D_rho*R_dRho*D_rho';
        z_ddPhi=ddPhi-hat_ddPhi;
%         outlier=abs(z_ddPhi)>1.5*mean(z_ddPhi);
%         z_ddPhi(outlier)=[];
         D_phi=D;
%         R_ddPhi=D_phi*R_dPhi*D_phi';
        %% Weighted-Least-Squared-Error Estimation
        u=data(:,9:11);
        pref=hatp;
        Nii=size(data,1);
        H=[-D_phi*u zeros(numel(z_ddPhi),3) lambda1*D_phi;...
            -D_rho*u zeros(numel(z_ddRho),Nii+3)];
        z=[z_ddPhi;z_ddRho];
        R=blkdiag(R_ddPhi,R_ddRho);
        hatN=data(:,15);
        hatx=[hatp;hatv;hatN]; %ith-prior
        x0=[hatx0;zeros(Nii,1)]; %first-prior
        Pna0=zeros(Nii);
        %% Cycle-slip
        cycleSlip(k)=sum(data(:,12));
        %%
        for s=1:Nii
            %        res_rho(sats==data(s,1),k)=z_ddRho(s);
            if s<Nii
                res_phi(sats==data(s,1),k)=z_ddPhi(s);
            end
            satsOnView(sats==data(s,1),k)=1;
            if k==1||(k>1 && satsOnView(sats==data(s,1),k-1)==0)
                satsTime(sats==data(s,1))=time(k);
                Pna(sats==data(s,1),sats==data(s,1),k)=pN0;
            end
            
            if data(s,12)~=0
                %cycle-slip
                Pna(sats==data(s,1),sats==data(s,1),k)=pN0;
            end
            
            Pna0(s,s)=Pna(sats==data(s,1),sats==data(s,1),k);
            x0(6+s)=Na0(sats==data(s,1));
        end
       
        %%
        P0=blkdiag(Px(:,:,k),Pna0);
        K=P0*H'/(R+H*P0*H');
        hatx=hatx+K*z;
        hatp=hatx(1:3);
        hatv=hatx(4:6);
        err=norm(hatp-pref);
        iter=iter+1;
    end
    %% result
    hx(:,k)=hatx(1:6);
    P=(eye(size(P0,1))-K*H)*P0;
    if min(eig(P))<0
        warning('Covariance non-positive!')
    end
    
    Px(:,:,k)=P(1:6,1:6);
    for s=1:Nii
        hatNa(sats==data(s,1),k)=round(hatx(6+s));
        Pna(sats==data(s,1),sats==data(s,1),k)=P(6+s,6+s);
    end
    %% Ambiguity resolution
    hatN=D*hatx(7:end);
    G=blkdiag(eye(6),D);
    W=G*P*G';
    WN=W(7:end,7:end);
    [optN,rr] = mlambda(WN,hatN,2);
    Ra(k)=rr(2)/rr(1);
    Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
    if Ra(k)>=Rthres && sqrt(trace(Pn))<0.5 && flagFix==1
        %% Update Ambiguities for the EKF
        Nii=size(hatN,1);
        GN=[zeros(Nii,6) D];
        K=P*GN'/(sigmaN^2*eye(Nii)+GN*P*GN');
        hatx=hatx+K*(optN(:,1)-GN*hatx);
        P=P-K*(sigmaN^2*eye(Nii)+GN*P*GN')*K';
        hx(:,k)=hatx(1:6);
        for s=1:Nii
            hatNa(sats==data(s,1),k)=hatx(6+s);
            Pna(sats==data(s,1),sats==data(s,1),k)=P(6+s,6+s);
            sigmaNp(sats==data(s,1))=sigmaNfix;
        end
    else
        for s=1:Nii
            sigmaNp(sats==data(s,1))=sigmaNfloat;
        end
    end
    %% Compute the Dilusion of Precision (DOP)
    Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
    dop(k,1)=sqrt(trace(Pn)); %Geometric (GDOP)
    dop(k,2)=sqrt(Pn(1,1)+Pn(2,2)); %Horizontal (HDOP)
    dop(k,3)=sqrt(Pn(3,3)); %Vertical (VDOP)
    trPna(k)=sqrt(trace(Pna(:,:,k))); %Ambiguity trace
    if ~mod(k,round(length(time)/10))
        fprintf('*')
    end
    %% time-update (prediction)
    if k<N
        [hx(:,k+1),Px(:,:,k+1)]=filterPred(hx(:,k),Px(:,:,k),dt,Cen*Qvn*Cen');
        hatNa(:,k+1)=hatNa(:,k);
%         hx(:,k+1)=hx(:,k);
%         Px(:,:,k+1)=Px(:,:,k);
        Pna(:,:,k+1)=Pna(:,:,k);%+10*diag(sigmaNp);
    end
end
%% solution
hx=hx';
fprintf('\nKinematic processing finished!\n')
% fprintf('Max. baseline: %.2fkm\n',max(baseline)/1e3)
% if max(baseline)/1e3 >10
%     warning('Long-baseline!')
% end
res.rho=res_rho;
res.phi=res_phi;
% hatNa(hatNa==0)=NaN;
% res.rho(res.rho==0)=NaN;
% res.phi(res.phi==0)=NaN;
% satsOnView(satsOnView==0)=NaN;
%% output
out.time=time;
out.hx=hx;
out.hatNa=hatNa;
out.satsOnView=satsOnView;
out.sats=sats;
out.Ra=Ra;
out.cycleSlip=cycleSlip;
out.res=res;
out.refSatId=refSatId;
out.Px=Px;
out.Pna=Pna;
end