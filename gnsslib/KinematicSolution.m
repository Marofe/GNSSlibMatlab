function [time,hx,hatNa,satsOnView,sats,Ra,cycleSlip]=KinematicSolution(ephemeris,obsBase,obsRover,pbase,p0rover,elevMask)
% Generate the high precision solution (~5cm) for GPS constellation
% input:
%       nav         -> broadcast nav messages
%       obs         -> observables (gps only)
%       p0          -> approximated initial position (ECEF)
%       elevMask    -> minimal elevation angle
fprintf('\nProcessing Kinematic Mode...\n')
%% constants
c=2.99792458e8; %Speed of light in vacuum (m/s)
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
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
Pna=zeros(Nsats,N); %single-difference ambiguity variance
% hxs=zeros(6,N);
% Pxs=zeros(6,6,N);
% hatNas=zeros(Nsats,N); %single-difference ambiguity
% Pnas=zeros(Nsats,N); %single-difference ambiguity variance
dop=zeros(N,3);
cycleSlip=zeros(N,1);
satsOnView=zeros(Nsats,N); %sats on
%% Params
iterMax=10; %maximum iterations of the ILS-update step
dt=1; %1hz GPS
sigmaVn=1; %north-acc noise (m/s/sqrt(s))
sigmaVe=1; %east-acc noise (m/s/sqrt(s))
sigmaVd=2; %down-acc noise (m/s/sqrt(s))
sigmaN=0.0001; %Ambiguity fix
sigmaNfloat=1;
sigmaNfix=1e-3;
sigmaNp=sigmaNfloat*ones(Nsats,1); %Ambiguity fix
Rthres=3;
%code/phase std = a*SNR+b/sin(elev)+c
%code
a1=0.1; %SNR contribution 
b1=.01; %Elev contribution 
c1=.5; %Independent noise 
%phase
a2=0; %SNR contribution 
b2=0.003; %Elev contribution 
c2=0.003; %Independent noise 
%Initial variance
p0=[1000*ones(1,3),10*ones(1,3)]; %position & velocity initial variance
pN0=10000; %initial ambiguity & cycle slip variance
flagFix=1;
%% Dynamic model (Const. Velocity)
A=[eye(3) dt*eye(3);zeros(3) eye(3)];
Qv=diag([sigmaVn^2,sigmaVe^2,sigmaVd^2])*dt;
%% Initialization
hx(:,1)=[p0rover;zeros(3,1)];
hatNa(:,1)=zeros(Nsats,1);
Px(:,:,1)=diag(p0);
Pna(:,1)=pN0*ones(Nsats,1);
%% Iterated-EKF-Processing
for k=1:N
    iter=1;
    err=inf;
    hatp=hx(1:3,k); %prior position
    hatv=hx(4:6,k); %prior velocity
    hatx0=hx(:,k);
    Px0=Px(:,:,k);
    %% Epochs
    %time(k) => time of epoch (toe)
    epoch=obsRover(obsRover(:,2)==time(k),:);
    epochBase=obsBase(obsBase(:,2)==time(k),:);
    % obs=[gpsWeek,tow,satConst,satID,pseudorange,rangeLLI,rangeStrength,phase,phaseLLI,phaseStrength,doppler,SNR];
    satInfo=getSatellitesInfo(epoch,epochBase,nav,hatp);
    % satInfo=[satID flagMissSat elev azim psat' L']
    %% take as reference sat the highest elevation
    [maxElev,refSat]=max(satInfo(satInfo(:,2)==0,3));
    refSatId=satInfo(refSat,1);
    %satref position
    psatRef=satInfo(refSat,5:7)';
    u1=satInfo(refSat,8:10);
    satsOnView(sats==refSatId,k)=1;
    %base obs sat ref
    baseRef=epochBase(epochBase(:,4)==refSatId,:);
    roverRef=epoch(refSat,:);
    refSatInfo=satInfo(refSat,:);
    %remove obs of refsat from epoch and satInfo
    epoch(refSat,:)=[];
    satInfo(refSat,:)=[];
    %%
    %rover ref measurements
    rho1_r=roverRef(5); %pseudo-range sat1 (rover)
    Phi1_r=roverRef(8)*lambda1; %phase-range sat1 (rover)
    SNR1_r=roverRef(end); %SNR sat1 (rover)
    
    sigmaR1_rho=a1*SNR1_r+b1/sin(maxElev)+c1;
    sigmaR1_phi=a2*SNR1_r+b2/sin(maxElev)+c2;
    
    %base ref measurements
    rho1_b=baseRef(5); %pseudo-range refsat (base)
    Phi1_b=baseRef(8)*lambda1; %phase-range refsat (base)
    SNR1_b=baseRef(end); %SNR refsat (base)
    
    %Assumption: elevBase=elevRover
    sigmaR1_rho_b=a1*SNR1_b+b1/sin(maxElev)+c1;
    sigmaR1_phi_b=a2*SNR1_b+b2/sin(maxElev)+c2;
    
    R1_rho=sigmaR1_rho^2+sigmaR1_rho_b^2;
    R1_phi=sigmaR1_phi^2+sigmaR1_phi_b^2;
    
    Ni=size(epoch,1); %number of satellites (=measurements)
    %% Prior Ambiguities
    if k==1
        Na0=zeros(Nsats,1);
    else
        hatNa(sats==refSatId,k)=hatNa(sats==refSatId,k-1);
        Pna(sats==refSatId,k)=Pna(sats==refSatId,k-1);
        for s=1:Ni
            hatNa(sats==epoch(s,4),k)=hatNa(sats==epoch(s,4),k-1);
            Pna(sats==epoch(s,4),k)=Pna(sats==epoch(s,4),k-1);
        end
        Na0=hatNa(:,k);
    end
    Pna0=Pna(:,k);
    %% ILS
    while err>1e-3 && iter<=iterMax
        lla=SingleLlaFromEcef(hatp);
        Cen=DCM_en(lla(1),lla(2));
        z1=[];
        z2=[];
        satsOn=refSatId;
        u=u1;
        R_rho=R1_rho; %single-diff
        R_phi=R1_phi; %single-diff
        
        for m=1:Ni
            elev=satInfo(m,3);
            if rad2deg(elev)>=elevMask
                satId=epoch(m,4);
                base=epochBase(epochBase(:,4)==satId,:);
                if size(base,1)>0 %if satId was found in the base epoch
                    %% Cycle-slip
                    if epoch(m,9)~=0
                        cycleSlip(k)=cycleSlip(k)+1;
                        Pna0(sats==epoch(m,4))=pN0;
                        hatNa(sats==epoch(m,4))=0;
                        Na0(sats==epoch(m,4))=0;
                        sigmaNp(sats==epoch(m,4))=sigmaNfloat;
                    end
                    %% Observables
                    %rover
                    rho_r=epoch(m,5); %pseudo-range 
                    Phi_r=epoch(m,8)*lambda1; %phase-range 
                    %doppler=-epoch(m,11)*lambda1; %range-rate
                    SNR_r=epoch(m,end); %Signal-to-Noise ratio
                    
                    sigmaR_rho=a1*SNR_r+b1/sin(elev)+c1;
                    sigmaR_phi=a2*SNR_r+b2/sin(elev)+c2;
%                     if epoch(m,9)==2
%                         sigmaR_phi=1000*sigmaR_phi;
%                     end
                    
                    %base
                    rho_b=base(5); %pseudo-range
                    Phi_b=base(8)*lambda1; %phase-range
                    %doppler=-epoch(m,11)*lambda1; %range-rate
                    SNR_b=base(end);  %Signal-to-Noise ratio
                    
                    sigmaR_rho_b=a1*SNR_b+b1/sin(elev)+c1;
                    sigmaR_phi_b=a2*SNR_b+b2/sin(elev)+c2;
                    
%                     if base(9)==2
%                         sigmaR_phi_b=1000*sigmaR_phi_b;
%                     end
                    
                    Rk_rho=sigmaR_rho^2+sigmaR_rho_b^2; %(single-difference)
                    Rk_phi=sigmaR_phi^2+sigmaR_phi_b^2; %(single-difference)
                    
                    %% Compute satellite position (at emission time)
                    psat=satInfo(m,5:7)';
                    range=psat-hatp;
                    bar_rho=norm(range); %approx. geometric range
                    %% Compute Double Difference
                    drho1=rho1_r-rho1_b;
                    drho=rho_r-rho_b;
                    ddrho=drho1-drho; %rho_rb_1j
                    
                    dphi1=Phi1_r-Phi1_b;
                    dphi=Phi_r-Phi_b;
                    ddPhi=dphi1-dphi; %rho_rb_1j
                    %% Pseudo-range
                    bar_drho1=norm(psatRef-hatp)-norm(psatRef-pbase);
                    bar_drho=norm(psat-hatp)-norm(psat-pbase);
                    bar_ddrho=bar_drho1-bar_drho;
                    
                    hat_ddrho=bar_ddrho; %predicted pseudo-range
                    hat_ddPhi=bar_ddrho+lambda1*(hatNa(sats==refSatId,k)-hatNa(sats==satId,k));  %predicted phase-range
                    %% Residuals
                    z_drho=drho-bar_drho;
                    z_dphi=dphi-bar_drho-lambda1*hatNa(sats==satId,k);
                    z1=[z1;ddPhi-hat_ddPhi]; %phase-range residual
                    z2=[z2;ddrho-hat_ddrho]; %pseudo-range residual
                    u=[u;range'/bar_rho];
                    satsOnView(sats==satId,k)=1;
                    satsOn=[satsOn;satId];
                    R_rho=blkdiag(R_rho,Rk_rho);
                    R_phi=blkdiag(R_phi,Rk_phi);
                else
                    warning('satelite not found in the base obs!')
                end
            end
        end
        %% Weighted-Least-Squared-Error Estimation
        pref=hatp;
        Nii=sum(satsOnView(:,k)==1)-1;
        D=[ones(Nii,1) -eye(Nii)];
        H=[-D*u zeros(Nii,3) lambda1*D;...
           -D*u zeros(Nii,Nii+1+3)];
        z=[z1;z2];
        R=blkdiag(D*R_phi*D',D*R_rho*D');
        hatx=[hatp;hatv;zeros(Nii+1,1)]; %ith-prior
        x0=[hatx0;zeros(Nii+1,1)]; %first-prior
        Pnai0=zeros(Nii+1,1);
        
        for s=1:Nii+1
            hatx(6+s)=hatNa(sats==satsOn(s),k);
            x0(6+s)=Na0(sats==satsOn(s));
            Pnai0(s)=Pna0(sats==satsOn(s));
            if Pnai0(s)==0
                Pnai0(s)=pN0;
            end
            Pnai0(s)=Pnai0(s)+sigmaNp(sats==satsOn(s));
        end
        
        P0=blkdiag(Px0,diag(Pnai0));
        K=P0*H'/(R+H*P0*H');
        hatx=x0+K*(z+H*(hatx-x0));
        hatp=hatx(1:3);
        hatv=hatx(4:6);
        err=norm(hatp-pref);
        if err<1e3
            for s=1:Nii+1
                hatNa(sats==satsOn(s),k)=hatx(6+s);
            end
        end
        iter=iter+1;
    end
    %% result
    hx(:,k)=hatx(1:6);
    P=P0-K*(R+H*P0*H')*K';
    Px(:,:,k)=P(1:6,1:6);
    for s=1:Nii+1
        Pna(sats==satsOn(s),k)=P(6+s,6+s);
    end
    %% Ambiguity resolution
    hatN=D*hatx(7:end);
    G=blkdiag(eye(6),D);
    W=G*P*G';
    WN=W(7:end,7:end);
%     iWN=eye(size(WN,1))/WN;
%     WXN=W(1:6,7:end);
    [optN,rr] = mlambda(WN,hatN,2);
    Ra(k)=rr(2)/rr(1);
    Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
    if Ra(k)>=Rthres && sqrt(trace(Pn))<0.5 && flagFix==1
%         hxx=hx(:,k)+WXN*iWN*(optN(:,1)-hatN);
        %% Update Ambiguities for the EKF
        Nii=size(hatN,1);
        GN=[zeros(Nii,6) D];
        K=P*GN'/(sigmaN^2*eye(Nii)+GN*P*GN');
        hatx=hatx+K*(optN(:,1)-GN*hatx);
        %hatx(7:end)=round(hatx(7:end));
        P=P-K*(sigmaN^2*eye(Nii)+GN*P*GN')*K';
        hx(:,k)=hatx(1:6);
        Pnai=P(7:end,7:end);      
        Pnai=diag(Pnai);
        for s=1:Nii
            hatNa(sats==satsOn(s),k)=hatx(6+s);
            Pna(sats==satsOn(s),k)=Pnai(s);
            sigmaNp(sats==satsOn(s))=sigmaNfix;
        end
    else
        for s=1:Nii
            sigmaNp(sats==satsOn(s))=sigmaNfloat;
        end
    end 
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
        hx(:,k+1)=A*hx(:,k);
        Q=blkdiag(zeros(3),Cen*Qv*Cen');
        Px(:,:,k+1)=A*Px(:,:,k)*A'+Q;
    end
end
%% Smoothing
% hatxs=[hx(:,end);hatNa(:,end)];
% hxs(:,end)=hx(:,end);
% Ps=blkdiag(Px(:,:,end),diag(Pna(:,end)));
% for k=N-1:-1:1
%     hatx0=[hx(:,k);hatNa(:,k)];
%     P0=blkdiag(Px(:,:,k),diag(Pna(:,k)));
%     F=blkdiag(A,eye(Nsats));
%     hatx=F*hatx0;
%     P=F*P0*F'+blkdiag(Q,5e-3*eye(Nsats));
%     G=P0*F'/P;
%     hatxs=hatx0+G*(hatxs-hatx);
%     Ps=P0+G*(Ps-P)*G';
%     hxs(1:6,k)=hatxs(1:6);
%     Pxs(:,:,k)=Ps(1:6,1:6);
%     hatNas(:,k)=hatxs(7:end);
%     Pnas(:,k)=diag(Ps(7:end,7:end));
% end
%% solution
hx=hx';
fprintf('\nKinematic processing finished!\n')
% fprintf('Max. baseline: %.2fkm\n',max(baseline)/1e3)
% if max(baseline)/1e3 >10
%     warning('Long-baseline!')
% end
end