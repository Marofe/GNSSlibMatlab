function [time,hx,satsOnView,sats,Ra]=KinematicSolution2(ephemeris,obsBase,obsRover,pbase,p0rover,elevMask)
% Generate the high precision solution (~5cm) for GPS constellation
% input:
%       nav         -> broadcast nav messages
%       obs         -> observables (gps only)
%       p0          -> approximated initial position (ECEF)
%       elevMask    -> minimal elevation angle
fprintf('\nProcessing Kinematic Mode...\n')
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
hx=zeros(6+Nsats,N);
Px=zeros(6+Nsats,6+Nsats,N);
satsOnView=zeros(Nsats,N); %sats on
Ra=zeros(N,1); %AR-Ratio
%% Params
iterMax=10; %maximum iterations of the ILS-update step
dt=1; %1hz GPS
sigmaVn=1; %north-acc noise (m/s/sqrt(s))
sigmaVe=1; %east-acc noise (m/s/sqrt(s))
sigmaVd=2; %down-acc noise (m/s/sqrt(s))
sigmaN=0.0001; %Ambiguity fix-hold
sigmaNp=sqrt(1); %Ambiguity 
Rthres=3; %ratio-threshould
trPxthres=4; %trace Px threshould
flagSmooth=0;
%code/phase std = a*SNR+b/sin(elev)+c
%code
a1=0.1; %SNR contribution 
b1=1; %Elev contribution 
c1=1; %Independent noise 
%phase
a2=0; %SNR contribution 
b2=0.003; %Elev contribution 
c2=0.003; %Independent noise 
%Initial variance
p0=1000; %position & velocity initial variance
pN0=10000; %initial ambiguity & cycle slip variance
%% Dynamic model (Const. Velocity)
A=[eye(3) dt*eye(3);zeros(3) eye(3)];
F=blkdiag(A,eye(Nsats));
Qv=diag([sigmaVn^2,sigmaVe^2,sigmaVd^2]);
%% Initialization
hx(:,1)=[p0rover;zeros(3+Nsats,1)];
Px(:,:,1)=blkdiag(p0*eye(6),pN0*eye(Nsats)); %initial covariance (position&&velocity);
%% Iterated-EKF-Processing
for k=1:N
    iter=1;
    err=inf;
    hatp=hx(1:3,k); %prior position
    hatv=hx(4:6,k); %prior velocity
    hatN=hx(7:end,k); %prior ambiguity
    hatx0=hx(:,k);
    P0=Px(:,:,k);
    %% Epochs
    %time(k) => time of epoch (toe)
    %select all measurements for t=time(k)
    epochRover=obsRover(obsRover(:,2)==time(k),:); 
    epochBase=obsBase(obsBase(:,2)==time(k),:);
    % obs=[gpsWeek,tow,satConst,satID,pseudorange,rangeLLI,rangeStrength,phase,phaseLLI,phaseStrength,doppler,SNR];
    % Compute all sat position within the k-th epoch
    satInfo=getSatellitesInfo(epochRover,epochBase,nav,hatp);
    %satInfo=[satId flagMissSat elev azim psat' L'];
    %% take as reference the highest elevation sat 
    [maxElev,refSat]=max(satInfo(satInfo(:,2)==0,3));
    refSatId=satInfo(refSat,1);
    refSatInfo=satInfo(refSat,:);
    %satref position
    psatRef=satInfo(refSat,5:7)';
    u1=satInfo(refSat,8:10);
    satsOnView(sats==refSatId,k)=1;
    %base obs sat ref
    baseRef=epochBase(epochBase(:,4)==refSatId,:);
    roverRef=epochRover(refSat,:);
    %remove obs of refsat from epoch and satInfo
    epochRover(refSat,:)=[];
    satInfo(refSat,:)=[];
    %% Ref measurements
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
    
    Ni=size(epochRover,1); %number of satellites (=measurements)  
    %% ILS
    while err>1e-3 && iter<=iterMax
        lla=SingleLlaFromEcef(hatp);
        Cen=DCM_en(lla(1),lla(2));
        z1=[];
        z2=[];
        satsOn=[];
        R_phi=R1_phi; %covariance of the phase-range (DD)
        R_rho=R1_rho; %covariance of the code-range (DD)
        u=[];
        for m=1:Ni
            elev=satInfo(m,3);
            if rad2deg(elev)>=elevMask
                satId=epochRover(m,4);
                base=epochBase(epochBase(:,4)==satId,:);
                if size(base,1)>0  %if satId was found in the base epoch
                    %% Observables
                    %rover
                    rho_r=epochRover(m,5); %pseudo-range
                    Phi_r=epochRover(m,8)*lambda1; %phase-range
                    %doppler=-epoch(m,11)*lambda1; %range-rate
                    SNR_r=epochRover(m,end); %Signal-to-Noise ratio
                    
                    sigmaR_rho=a1*SNR_r+b1/sin(elev)+c1;
                    sigmaR_phi=a2*SNR_r+b2/sin(elev)+c2;
                    
                    %base
                    rho_b=base(5); %pseudo-range (~5m) (base)
                    Phi_b=base(8)*lambda1; %phase-range (~10cm+ambiguity) (base)
                    %doppler=-epoch(m,11)*lambda1; %range-rate
                    SNR_b=base(end);  %Signal-to-Noise ratio (base)
                    
                    sigmaR_rho_b=a1*SNR_b+b1/sin(elev)+c1;
                    sigmaR_phi_b=a2*SNR_b+b2/sin(elev)+c2;
                    
                    %% Side-slip
                    if epochRover(m,9)==2 || epochBase(m,9)==2
                        id=sats==satId;
                        hatN(id)=0;
                        P0(id,id)=pN0;
                    end
                    
                    %% Satellite position (at emission time)
                    psat=satInfo(m,5:7)';
                    range=psat-hatp;
                    bar_rho=norm(range); %approx. geometric range
                    %% Compute Double Difference
                    ddrho=(rho1_r-rho1_b)-(rho_r-rho_b); %rho_rb_1j
                    ddPhi=(Phi1_r-Phi1_b)-(Phi_r-Phi_b); %rho_rb_1j
                    %% Pseudo-range
                    bar_ddrho=(norm(psatRef-hatp)-norm(psatRef-pbase))-(norm(psat-hatp)-norm(psat-pbase));
                    hat_ddrho=bar_ddrho; %predicted pseudo-range
                    hat_ddPhi=bar_ddrho+lambda1*(hatN(sats==refSatId)-hatN(sats==satId));  %predicted phase-range
                    %% Residuals
                    z_phi=ddPhi-hat_ddPhi;
                    z_rho=ddrho-hat_ddrho;
                    if z_rho/sigmaR_rho<=3
                    z1=[z1;z_phi]; %phase-range residual
                    z2=[z2;z_rho]; %pseudo-range residual
                    satsOnView(sats==satId,k)=1;
                    u=[u;u1-satInfo(m,8:10)]; %r1-rj=-(u1-uj)
                    satsOn=[satsOn;satId];
                    Rk_rho=sigmaR_rho^2+sigmaR_rho_b^2; %(single-difference)
                    Rk_phi=sigmaR_phi^2+sigmaR_phi_b^2; %(single-difference)
                    R_phi=[R_phi;Rk_phi];
                    R_rho=[R_rho;Rk_rho];
                    else
                        warning('outlier!')
                    end
                else
                    warning('satelite not found in the base obs!')
                end
            end
        end
        %% Weighted-Least-Squared-Error Estimation
        pref=hatp;
        Nii=sum(satsOnView(:,k)==1)-1;
        if Nii~=Ni
            warning('Error in the sat counting... Nii~Ni!')
        end
        D=zeros(Nii,Nsats);
        D(:,sats==refSatId)=1;
        for jj=1:Nii
            D(jj,sats==satsOn(jj))=-1;
        end
        H=[-u zeros(Nii,3) lambda1*D;-u zeros(Nii,3+Nsats)];
        z=[z1;z2];
        
        R_rho=diag(R_rho);
        R_phi=diag(R_phi);
        D0=[ones(Nii,1) -eye(Nii)];
        
        R=blkdiag(D0*R_phi*D0',D0*R_rho*D0');
        hatx=[hatp;hatv;hatN]; %ith-prior
        %hatx0%first-prior
        K=P0*H'/(R+H*P0*H');
        hatx=hatx0+K*(z+H*(hatx-hatx0));
        err=norm(hatx(1:3)-pref);
        if err<1e3
            %inlier
            hatp=hatx(1:3);
            hatv=hatx(4:6);
            hatN=hatx(7:end);
        end
        iter=iter+1;
    end
    %% result
    hx(:,k)=hatx;
    P=P0-K*(R+H*P0*H')*K';
    Px(:,:,k)=P;
    trP(k)=trace(P);
    trPx(k)=trace(P(1:3,1:3));
    trPv(k)=trace(P(4:6,4:6));
    trNa(k)=trace(P(7:end,7:end));
    
    %% Ambiguity resolution
    ddN=D*hatN;
    G=blkdiag(eye(6),D);
    W=G*P*G';
    WN=W(7:end,7:end);
    iWN=eye(size(WN,1))/WN;
    WXN=W(1:6,7:end);
    [optN,rr] = mlambda(WN,ddN,2);
    Ra(k)=((ddN-optN(:,2))'*iWN*(ddN-optN(:,2)))/((ddN-optN(:,1))'*iWN*(ddN-optN(:,1)));
    if Ra(k)>=Rthres && trPx(k)<trPxthres
        %% Update Ambiguities for the EKF
        Nii=size(ddN,1);
        GN=[zeros(Nii,6) D];
        K=P*GN'/(sigmaN^2*eye(Nii)+GN*P*GN');
        hatx=hatx+K*(optN(:,1)-GN*hatx);
        hatx(7:end)=round(hatx(7:end));
%         P=P-K*(sigmaN^2*eye(Nii)+GN*P*GN')*K';
        hx(:,k)=hatx;
    end
    
    %% Compute the Dilusion of Precision (DOP)
%     Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
%     dop(k,1)=sqrt(trace(Pn)); %Geometric (GDOP)
%     dop(k,2)=sqrt(Pn(1,1)+Pn(2,2)); %Horizontal (HDOP)
%     dop(k,3)=sqrt(Pn(3,3)); %Vertical (VDOP)
    if ~mod(k,round(length(time)/10))
        fprintf('*')
    end
    
    %% time-update (prediction)
    if k<N
        hx(:,k+1)=F*hx(:,k);
        Q=blkdiag(zeros(3),Cen*Qv*Cen'*dt,sigmaNp^2*eye(Nsats));
        Px(:,:,k+1)=F*Px(:,:,k)*F'+Q;
        Ra(k+1)=Ra(k);
    end
end
%% Smoothing
if flagSmooth==1
hatxs=hx(:,end);
Ps=Px(:,:,end);
fprintf('\nSmoothing...')
for k=N-1:-1:1
    hatx0=hx(:,k);
    P0=Px(:,:,k);
    hatx=F*hatx0;
    lla=SingleLlaFromEcef(hatx(1:3));
    Cen=DCM_en(lla(1),lla(2));
    Q=blkdiag(zeros(3),Cen*Qv*Cen'*dt,sigmaNp^2*eye(Nsats));
    P=F*P0*F'+Q;
    G=P0*F'/P;
    hatxs=hatx0+G*(hatxs-hatx);
    Ps=P0+G*(Ps-P)*G';
    hx(:,k)=hatxs;
    Px(:,:,k)=Ps;
    %% Ambiguity
    dN=hatxs(7:end);
    WN=Ps(7:end,7:end);
    iWN=eye(size(WN,1))/WN;
    WXN=Ps(1:6,7:end);
    [optN,rr] = mlambda(WN,dN,2);
    Ra(k)=((dN-optN(:,2))'*iWN*(dN-optN(:,2)))/((dN-optN(:,1))'*iWN*(dN-optN(:,1)));
    if Ra(k)>=Rthres
        %% Update Ambiguities for the EKF
        GN=[zeros(Nsats,5) eye(Nsats)];
        K=P*GN'/(sigmaN^2*eye(Nsats)+GN*P*GN');
        hatxs=hatxs+K*(optN(:,1)-GN*hatxs);
        hatxs(7:end)=round(hatxs(7:end));
        Ps=Ps-K*(sigmaN^2*eye(Nsats)+GN*Ps*GN')*K';
        hx(:,k)=hatxs;
        Px(:,:,k)=Ps;
    end
    %%
    if ~mod(k,round(length(time)/10))
        fprintf('*')
    end
end
end
%% solution
% hx=hx';
fprintf('\nKinematic processing finished!\n')
% fprintf('Max. baseline: %.2fkm\n',max(baseline)/1e3)
% if max(baseline)/1e3 >10
%     warning('Long-baseline!')
% end
end