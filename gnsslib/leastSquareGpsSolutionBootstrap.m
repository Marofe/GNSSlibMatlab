function [time,p,ns,dop,atmDelay]=leastSquareGpsSolutionBootstrap(nav,obs,p0,elevMask,atmParam)
% Generate the standard precision solution (~5m) for GPS constellation
% input:
%       nav         -> broadcast nav messages
%       obs         -> observables (gps only)
%       p0          -> approximated initial position (ECEF)
%       elevMask    -> minimal elevation angle
%% Constants
mu=3986004.418e8;% gravitational constant (m^3/s^2)
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
R0=6378137; %Earth Radius
earth_ecc=0.0818191908425; %Eccentricity of Earth
phi_p=deg2rad(78.3); %magnet pole lat
lambda_p=deg2rad(291); %magnet pole lon
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
p=zeros(length(time),4);
dop=zeros(length(time),4);
ns=zeros(length(time),1);
atmDelay=zeros(length(time),2);
for k=1:length(time)
    %% Epoch
    trcv=time(k);
    epoch=obs(obs(:,2)==trcv,:);
    Ni=size(epoch,1);
    cdt=0;
    if k>1
        hatp=p(k-1,:)';
    else
        hatp=[p0;cdt];
    end
    s=1;
    pref=zeros(4,1);
    while norm(hatp-pref)>1e-3 && s<=15
        lla=SingleLlaFromEcef(hatp(1:3));
        res=[];
        H=[];
        W=[];
        totalTrop=0;
        totalIon=0;
        for m=1:Ni
            satId=epoch(m,4);
            id=find(nav(1,:)==satId);        %select all ephemeris for the satellite
            [~,j]=min(abs(trcv-nav(18,id)));  %seek for the most recent orbit parameters
            j=id(j);
            toe=nav(18,j); %time of ephemeris
            doy=doyFromGPST(epoch(m,1:2));
            TGD=nav(22,j); %Total Group Delay (instrumental error)
            %% Observables
            r=epoch(m,5); %pseudo-range (~5m)
            %phase=epoch(m,8); %carrier-phase (~10cm+ambiguity)
            %doppler=epoch(m,11); %doppler frequency
            SNR=epoch(m,end); %Signal-to-Noise ratio
            %% Compute satellite clock offset
            af0=nav(19,j);
            af1=nav(20,j);
            af2=nav(2,j);
            tc=toe;
            tems0=trcv+(cdt-r)/c;
            dtsat0=af0+af1*(tems0-tc)+af2*(tems0-tc)^2;
            %% Compute relativistic correction
            [~,E,a]=satPosition(nav(:,j),tems0);
            ecc=nav(6,j);
            dtrel0=-2*sqrt(mu*a)/(c^2)*ecc*sin(E);
            %% Compute the emission time
            dTs=dtsat0+dtrel0-TGD;
            tems=tems0-dTs; %(GPST)
            %update satellite clock offset
            dtsat=af0+af1*(tems-tc)+af2*(tems-tc)^2;
            %% Compute satellite position (emission time)
            % From Broadcast nav msg
            [satPos0,E,a]=satPosition(nav(:,j),tems);
            range0=satPos0-hatp(1:3);
            rho0=norm(range0); %approx. geometric range
            %% Compute satellite position (reception time)
            % Take in account the Earth's rotation
            theta=wie*rho0/c;
            satPos=rotZ(theta)*satPos0;
            rhoSat=norm(satPos);
            range=satPos-hatp(1:3);
            rho=norm(range); %approx. geometric range
            rhoRcv=norm(hatp(1:3)); %predicted rcv position
            Cen=DCM_en(lla(1),lla(2));
            los=Cen'*range/rho; %Line of Sight (LOS)
            azim=atan2(los(1),los(2)); %Azimute (NED frame)
            elev=asin(-los(3)); %Elevation (rad) (NED frame)
            if rad2deg(elev)>elevMask
                %% update relativistic corrections
                drel=2*mu/(c^2)*log((rhoSat+rhoRcv+rho)/(rhoSat+rhoRcv-rho));
                dtrel=-2*sqrt(mu*a)/(c^2)*ecc*sin(E);
                %% Compute Tropospheric delay
                trop=troposphericModel(lla,doy,elev);
                totalTrop=totalTrop+trop;
                %% Compute Ionospheric delay
                %Klobuchar Coeff from nav msg.
                alpha0=atmParam(1,1);
                alpha1=atmParam(1,2);
                alpha2=atmParam(1,3);
                alpha3=atmParam(1,4);
                beta0=atmParam(2,1);
                beta1=atmParam(2,2);
                beta2=atmParam(2,3);
                beta3=atmParam(2,4);
                DC=5e-9;
                % Calculate the Earth centered angle
                lat=deg2rad(lla(1));
                lon=deg2rad(lla(2));
                RE=R0/sqrt(1-earth_ecc^2*sin(lat)^2);
                RR=RE/(RE+350000); %atm height~350km
                psi=pi/2-elev-asin(RR*cos(azim));
                % Compute the latitude of the IPP (Iono Pierce Point)
                phi_i=asin(sin(lat)*cos(psi)+cos(lat)*sin(psi)*cos(azim));
                % Compute the longitude of the IPP
                lambda_i=lon+psi*sin(azim)/cos(phi_i);
                % Find the geomagnetic latitude of the IPP
                phi_m=asin(sin(phi_i)*sin(phi_p)+cos(phi_i)*cos(phi_p)*cos(lambda_i-lambda_p));
                % Find the local time at the IPP
                t=trcv+cdt/c;
                if t>=86400
                    t=t-86400; %24h=86400s
                elseif t<0
                    t=t+86400;
                end
                ti=43200*lambda_i/pi+t; %12h=43200s
                % Compute the amplitude of ionospheric delay
                Ai=alpha0+alpha1*phi_m/pi+alpha2*phi_m^2/pi^2+alpha3*phi_m^3/pi^3;
                if Ai<0
                    Ai=0;
                end
                % Compute the period of ionospheric delay
                Pi=beta0+beta1*phi_m/pi+beta2*phi_m^2/pi^3+beta3*phi_m^3/pi^3;
                if Pi<72000
                    Pi=7200;
                end
                % Compute the phase of ionospheric delay
                Xi=2*pi*(ti-50400)/Pi; %14h=50400s
                % Compute the slant factor (ionospheric mapping function)
                F=(1-(RR*cos(elev))^2)^(-1/2);
                % Compute the ionospheric time delay
                if abs(Xi)<pi/2
                    ion=c*F*(DC+Ai*cos(Xi)); %day
                else
                    ion=c*F*DC; %night
                end
                totalIon=totalIon+ion;
                %% Pseudo-range
                rho=rho+drel; %geom dist + rel offset (~3cm)
                dTs=dtsat+dtrel-TGD; %total sat clock offset (sat+rel+inst)
                tilR=rho+cdt-c*dTs+trop+ion; %predicted pseudo-range
                res=[res;r-tilR];
                H=[H;-range'/rho 1];
                sigmaR=10/SNR+3*exp(-2*elev/(pi/2));
%                 a=0.13;b=0.53;c=10; %from MOPS
%                 sigmaR=a+b*exp(-rad2deg(elev)/c);
                W=[W 1/sigmaR^2];
            end
            emst(m)=tems;
        end
        %% Weighted-Least-Squared-Error Estimation
        W=diag(W);
        pref=hatp;
        hatp=hatp+(H'*W*H)\H'*W*res;
        cdt=hatp(4);
        s=s+1;
    end
    satValid=abs(res)<2.5;
    if sum(~satValid)~=0
        W=diag(W);
        W=diag(W(satValid));
        H=H(satValid,:);
        hatp=hatp+(H'*W*H)\H'*W*res(satValid);
    end
    %% Bootstrap (parametric)
    if k>1
        p00=p(k-1,:)';
    else
        p00=[p0; 0];
    end
    for w=1:1000
    hatpi(:,w)=epochEstimate(p00,epoch,nav,atmParam,time(k),hatp,emst,elevMask);
    end
    hatpm=mean(hatpi')';
    Pb=cov(hatpi');
    %% Result  
    p(k,:)=hatp';
    ns(k)=size(H,1);
    P=inv(H'*W*H);
    %% Compute the Dilusion of Precision (DOP)
    Pn=Cen'*P(1:3,1:3)*Cen; % ECEF->NED
    Pnb=Cen'*Pb(1:3,1:3)*Cen; % ECEF->NED
    dop(k,1)=sqrt(trace(Pn)); %Geometric (GDOP)
    dop(k,2)=sqrt(Pn(1,1)+Pn(2,2)); %Horizontal (HDOP)
    dop(k,3)=sqrt(Pn(3,3)); %Vertical (VDOP)
    dop(k,4)=sqrt(P(4,4)); %Time (TDOP)
    %
    dopb(k,1)=sqrt(trace(Pnb)); %Geometric (GDOP)
    dopb(k,2)=sqrt(Pnb(1,1)+Pnb(2,2)); %Horizontal (HDOP)
    dopb(k,3)=sqrt(Pnb(3,3)); %Vertical (VDOP)
    dopb(k,4)=sqrt(Pb(4,4)); %Time (TDOP)
    %% Atmospheric delay
    atmDelay(k,1)=totalTrop/ns(k); %ZTD (Zenith Tropospheric Delay)
    atmDelay(k,2)=totalIon/ns(k); %Ionospheic
    %
    if ~mod(k,round(length(time)/10))
        fprintf('*')
    end
end
fprintf('\nProcessing finished!')
end