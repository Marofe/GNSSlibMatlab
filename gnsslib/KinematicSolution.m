function [time,hx,hxs,hatNa,hatNas,baseline,satsOnView,sats]=KinematicSolution(ephemeris,obsBase,obsRover,pbase,p0rover,elevMask)
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
hxs=zeros(6,N);
Pxs=zeros(6,6,N);
hatNa=zeros(Nsats,N); %single-difference ambiguity
Pna=zeros(Nsats,N); %single-difference ambiguity variance
hatNas=zeros(Nsats,N); %single-difference ambiguity
Pnas=zeros(Nsats,N); %single-difference ambiguity variance
dop=zeros(N,3);
ns=zeros(N,1);
satsOnView=zeros(Nsats,N); %sats on
baseline=zeros(N,1);
%% Params
iterMax=15; %maximum iterations of the ILS-update step
dt=1; %1hz GPS
sigmaVn=0.1; %north-acc noise (m/s/sqrt(s))
sigmaVe=0.1; %east-acc noise (m/s/sqrt(s))
sigmaVd=0.1; %down-acc noise (m/s/sqrt(s))
Px0=10*eye(6); %initial covariance (position&&velocity)
Pna0=100*ones(Nsats,1); %initial covariance (ambiguities)
%% Dynamic model (Const. Velocity)
A=[eye(3) dt*eye(3);zeros(3) eye(3)];
Qv=diag([sigmaVn^2,sigmaVe^2,sigmaVd^2])*dt;
%% Initialization
hx(:,1)=[p0rover;zeros(3,1)];
hatNa(:,1)=zeros(Nsats,1);
Px(:,:,1)=Px0;
Pna(:,1)=Pna0;
%% Iterated-EKF-Processing
for k=1:N
    iter=1;
    err=inf;
    hatp=hx(1:3,k); %prior position
    hatv=hx(4:6,k); %prior velocity
    xx0=[hatp;hatv];
    %% Epochs
    %time(k) => time of epoch (toe)
    epoch=obsRover(obsRover(:,2)==time(k),:);
    epochBase=obsBase(obsBase(:,2)==time(k),:);
    % obs=[gpsWeek,tow,satConst,satID,pseudorange,rangeLLI,rangeStrength,phase,phaseLLI,phaseStrength,doppler,SNR];
    satInfo=getSatellitesInfo(epoch,epochBase,nav,hatp);
    %% take as reference sat the highest elevation
    [maxElev,refSat]=max(satInfo(satInfo(:,2)==0,3));
    refSatId=satInfo(refSat,1);
    %satref position
    psatRef=satInfo(refSat,5:7)';
    satsOnView(sats==refSatId,k)=1;
    %base obs sat ref
    baseRef=epochBase(epochBase(:,4)==refSatId,:);
    roverRef=epoch(refSat,:);
    %remove obs of refsat from epoch and satInfo
    epoch(refSat,:)=[];
    satInfo(refSat,:)=[];
    %%
    %rover ref measurements
    rho1_r=roverRef(5); %pseudo-range sat1 (rover)
    Phi1_r=roverRef(8)*lambda1; %phase-range sat1 (rover)
    SNR1_r=roverRef(end); %SNR sat1 (rover)
    
    sigmaR1=10/SNR1_r+3*exp(-2*maxElev/(pi/2));
    
    %base ref measurements
    rho1_b=baseRef(5); %pseudo-range refsat (base)
    Phi1_b=baseRef(8)*lambda1; %phase-range refsat (base)
    SNR1_b=baseRef(end); %SNR refsat (base)
    
    Ni=size(epoch,1); %number of satellites (=measurements)
    %% Prior Ambuities
    if k==1
        Na0=zeros(Nsats,1);
    else
        hatNa(sats==refSatId,k)=hatNa(sats==refSatId,k-1);
        for s=1:Ni
            hatNa(sats==satInfo(s,1),k)=hatNa(sats==satInfo(s,1),k-1);
            if epoch(s,9)~=0
                %side-slip
                Pna(sats==satInfo(s,1),k)=100;
            end
        end
        Na0=hatNa(:,k);
    end
    %% ILS
    while err>1e-3 && iter<=iterMax
        lla=SingleLlaFromEcef(hatp);
        Cen=DCM_en(lla(1),lla(2));
        z1=[];
        z2=[];
        satsOn=refSatId;
        u=satInfo(refSat,8:10);
        R=sigmaR1*sigmaR1';
        for m=1:Ni
            elev=satInfo(m,3);
            if rad2deg(elev)>=elevMask
                satId=epoch(m,4);
                base=epochBase(epochBase(:,4)==satId,:);
                if size(base,1)>0
                    %satId found in the base epoch
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
                    SNR_b=base(end);  %Signal-to-Noise ratio (base)
                    
                    %% Compute satellite position (at emission time)
                    psat=satInfo(m,5:7)';
                    range=psat-hatp;
                    bar_rho=norm(range); %approx. geometric range
                    %% Compute Double Difference
                    ddrho=(rho1_r-rho1_b)-(rho_r-rho_b); %rho_rb_1j
                    ddPhi=(Phi1_r-Phi1_b)-(Phi_r-Phi_b); %rho_rb_1j
                    %% Pseudo-range
                    bar_ddrho=(norm(psatRef-hatp)-norm(psatRef-pbase))-(norm(psat-hatp)-norm(psat-pbase));
                    hat_ddrho=bar_ddrho; %predicted pseudo-range
                    hat_ddPhi=bar_ddrho+lambda1*(hatNa(sats==refSatId,k)-hatNa(sats==satId,k));  %predicted phase-range
                    %% Residuals
                    z1=[z1;ddPhi-hat_ddPhi]; %phase-range residual
                    z2=[z2;ddrho-hat_ddrho]; %pseudo-range residual
                    u=[u;range'/bar_rho];
                    satsOnView(sats==satId,k)=1;
                    satsOn=[satsOn;satId];
                    sigmaR=10/SNR_r+3*exp(-2*elev/(pi/2));
                    R=blkdiag(R,sigmaR*sigmaR');
                else
                    warning('satelite not found in the base obs!')
                end
            end
        end
        %% Weighted-Least-Squared-Error Estimation
        pref=hatp;
        Nii=sum(satsOnView(:,k)==1)-1;
        D=[ones(Nii,1) -eye(Nii)];
        H=[-D*u zeros(Nii,3) lambda1*D;-D*u zeros(Nii,Nii+1+3)];
        z=[z1;z2];
        R=blkdiag(D*R/1000*D',D*R/10*D');
        hatx=[hatp;hatv;zeros(Nii+1,1)]; %ith-prior
        x0=[xx0;zeros(Nii+1,1)]; %first-prior
        Pnai=zeros(Nii+1,1);
        
        for s=1:Nii+1
            hatx(6+s)=hatNa(sats==satsOn(s),k);
            x0(6+s)=Na0(sats==satsOn(s));
            Pnai(s)=Pna(sats==satsOn(s),k);
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
                hatNa(sats==satsOn(s),k)=hatx(6+s);
            end
        end
        iter=iter+1;
    end
    %% result
    hx(:,k)=hatx(1:6);
    baseline(k)=norm(pbase-hx(1:3,k));
    P=P0-K*(R+H*P0*H')*K';
    Px(:,:,k)=P(1:6,1:6);
    Pnai=P(7:end,7:end);
    ns(k)=size(u,1);
    %% Ambiguity resolution
%         hatN=D*hatx(7:end);
%         G=blkdiag(eye(6),D);
%         W=G*P*G';
%         WN=W(7:end,7:end);
%         WXN=W(1:6,7:end);
%         optN=round(hatN);
%         hx(:,k)=hx(:,k)+WXN/WN*(optN-hatN);
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
        hx(:,k+1)=A*hx(:,k);
        Q=blkdiag(zeros(3),Cen*Qv*Cen');
        Px(:,:,k+1)=A*Px(:,:,k)*A'+Q;
        Pnai=Pnai+5e-3*eye(Nii+1);
        Pnai=diag(Pnai);
        for s=1:Nii+1
            Pna(sats==satsOn(s),k+1)=Pnai(s);
        end
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
fprintf('Max. baseline: %.2fkm\n',max(baseline)/1e3)
if max(baseline)/1e3 >10
    warning('Long-baseline!')
end
end