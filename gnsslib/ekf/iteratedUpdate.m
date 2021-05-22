function [hx,P,ns,totalIon,totalTrop,lla] = iteratedUpdate(hx0,P0,epoch,nav,atmParam,time,elevMask)
%% Constants
mu=3986004.418e8;% gravitational constant (m^3/s^2)
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
%%
Ni=size(epoch,1);
iter=1;
iterMax=15;
hx=hx0;
v=inf;
while norm(v)>1e-3 && iter<=iterMax
    lla=SingleLlaFromEcef(hx(1:3));
    z=[];
    H=[];
    W=[];
    totalIon=0;
    totalTrop=0;
    for m=1:Ni
        satId=epoch(m,4);
        id=find(nav(1,:)==satId);        %select all ephemeris for the satellite
        [~,j]=min(abs(time-nav(18,id)));  %seek for the most recent orbit parameters
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
        tems0=time+(hx(4)-r)/c;
        dtsat0=af0+af1*(tems0-tc)+af2*(tems0-tc)^2;
        %% Compute relativistic correction
        [~,E,a]=satPosition(nav(:,j),tems0);
        ecc=nav(6,j);
        dtrel0=-2*sqrt(mu*a)/(c^2)*ecc*sin(E);
        %% Compute the emission time
        dTs=dtsat0+dtrel0-TGD;
        tems=time+hx(4)/c-(r/c+dTs); %(GPST)
        %update satellite clock offset
        dtsat=af0+af1*(tems-tc)+af2*(tems-tc)^2;
        %% Compute satellite position (at emission time)
        % From Broadcast nav msg
        [satPos0,E,a]=satPosition(nav(:,j),tems);
        range0=satPos0-hx(1:3);
        rho0=norm(range0); %approx. geometric range
        %% Take in account the Earth's rotation
        theta=wie*rho0/c;
        satPos=rotZ(theta)*satPos0;
        rhoSat=norm(satPos);
        range=satPos-hx(1:3);
        rho=norm(range); %approx. geometric range
        rhoRcv=norm(hx(1:3)); %predicted rcv position
        Cen=DCM_en(lla(1),lla(2));
        los=Cen'*range/rho; %Line of Sight (LOS)
        azim=atan2(los(1),los(2)); %Azimute (NED frame)
        elev=asin(-los(3)); %Elevation (rad) (NED frame)
        elevd=rad2deg(elev); %Elevation (deg)
        if elevd>elevMask
            %% update relativistic corrections
            drel=2*mu/(c^2)*log((rhoSat+rhoRcv+rho)/(rhoSat+rhoRcv-rho));
            dtrel=-2*sqrt(mu*a)/(c^2)*ecc*sin(E);
            %% Compute Tropospheric delay
            trop=troposphericModel(lla,doy,elev);
            %% Compute Ionospheric delay
            %Klobuchar Coeff from nav msg.
            ion=ionosphericModel(atmParam,lla,elev,azim,time);
            totalIon=totalIon+ion;
            %% Pseudo-range
            rho=rho+drel; %geom dist + rel offset (~3cm)
            dTs=dtsat+dtrel-TGD; %total sat clock offset (sat+rel+inst)
            tilR=rho+hx(4)-c*dTs+trop+ion; %predicted pseudo-range
            z=[z;r-tilR];
            H=[H;-range'/rho 1];
            sigmaR=10/SNR+3*exp(-2*elev/(pi/2));
            W=[W sigmaR^2];
        end
    end
    %% EKF Update step (Iterated Version)
    R=diag(W);
    K=P0*H'/(R+H*P0*H');
    v=K*(z+H*(hx-hx0));
    hx=hx0+v;
    iter=iter+1;
end
%% Consistency test
% satValid=abs(z)<2.5;
% if sum(~satValid)~=0
%     W=diag(R);
%     R=diag(W(satValid));
%     H=H(satValid,:);
%     K=P0*H'/(R+H*P0*H');
%     v=K*(z(satValid)+H*(hx-hx0));
%     hx=hx0+v;
% end
%% result
hx=hx';
ns=size(H,1);
P=(eye(4)-K*H)*P0*(eye(4)-K*H)'+K*R*K'; %Joseph Form
end

