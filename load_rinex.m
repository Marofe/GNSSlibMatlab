close all
clear all
clc
format long
addpath('gnsslib')
%% Load RINEX V3.02 (Receiver INdependent EXchange File)
file='gnss_logger_dat-2021-01-27-13-14-20';
%file='gnss_logger_dat-2021-01-22-12-04-00';
ephemeris = readRinexNav([file '.nav']);
tic
[XYZ_station,allObs]=readRinex302([file '.obs']);
toc
%load('rinex_data20210122120400.mat')
%% Constants
mu=3986004.418e8;% gravitational constant (m^3/s^2)
we=7292115.0e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
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
%%
nav=ephemeris.gpsEphemeris;
satID=unique(nav(1,:));
gpsObs=allObs(:,3)==0;
obs=allObs(gpsObs,:);
time=unique(obs(:,2));
Ns=numel(satID);
for k=1:length(time)
    for m=1:Ns
        id=find(nav(1,:)==satID(m));        %select all ephemeris for the m-th satellite
        [~,j]=min(abs(time(k)-nav(18,id))); %seek the more recent orbit parameters for the m-th satellite
        j=id(j);                            %select ephemeris to compute orbit
        toe(k,m)=nav(18,j);
        M0(k,m)=nav(3,j);                        %mean anomaly
        sqrta(k,m)=nav(4,j);
        a=sqrta(k,m)^2;
        deltan(k,m)=nav(5,j);
        ecc(k,m)=nav(6,j);
        omega(k,m)=nav(7,j);
        cuc(k,m)=nav(8,j);
        cus(k,m)=nav(9,j);
        crc(k,m)=nav(10,j);
        crs(k,m)=nav(11,j);
        i0(k,m)=nav(12,j);
        doti0(k,m)=nav(13,j);
        cic(k,m)=nav(14,j);
        cis(k,m)=nav(15,j);
        Omega0(k,m)=nav(16,j);
        dotOmega0(k,m)=nav(17,j);
        %%
        t(k,m)=time(k)-toe(k,m);
        if t(k,m)>302400
            t(k,m)=t(k,m)-604800;
        elseif t(k,m)<-302400
            t(k,m)=t(k,m)+604800;
        end
        M=M0(k,m)+(sqrt(mu)/(sqrta(k,m)*a)+deltan(k,m))*t(k,m);
        E=M;
        while abs(E-M-ecc(k,m)*sin(E))>1e-15
            E=M+ecc(k,m)*sin(E);
        end
        %true anomaly
        v=atan(sqrt(1-ecc(k,m)^2)*sin(E)/(cos(E)-ecc(k,m)));
        %orbit latitude (true anomaly+arg perigee)
        u=omega(k,m)+v+cuc(k,m)*cos(2*(omega(k,m)+v))+cus(k,m)*sin(2*(omega(k,m)+v));
        r=a*(1-ecc(k,m)*cos(E))+crc(k,m)*cos(2*(omega(k,m)+v))+crs(k,m)*sin(2*(omega(k,m)+v));
        i=i0(k,m)+doti0(k,m)*t(k,m)+cic(k,m)*cos(2*(omega(k,m)+v))+cis(k,m)*sin(2*(omega(k,m)+v));
        Omega=Omega0(k,m)+(dotOmega0(k,m)-we)*t(k,m)-we*toe(k,m);
        C(:,:,k,m)=rotZ(-Omega)*rotX(-i)*rotZ(-u);
        satPos(:,k,m)=C(:,:,k,m)*[r;0;0];
    end
end
%%
gpsWeekStart = datetime(2021,1,27,'TimeZone','UTC');
plotEarth()
k=1;
%     timestamp = gpsWeekStart + seconds(time(k))
for m=1:Ns
plot3([0 satPos(1,k,m)],[0 satPos(2,k,m)],[0 satPos(3,k,m)],'r-','linewidth',1)
b=sqrt(1-ecc(k,m)^2)*sqrta(k,m)^2;
coord=ellipse3D(sqrta(k,m)^2,b,0,0,0,300,0,0,0);
coord=C(:,:,k,m)*coord;
plot3(coord(1,:),coord(2,:),coord(3,:),'w--')
plot3(satPos(1,k,m),satPos(2,k,m),satPos(3,k,m),'go','linewidth',2)
text(satPos(1,k,m),satPos(2,k,m),satPos(3,k,m),['G' num2str(satID(m))],'color','w')
satLegend{m}=['G' num2str(satID(m))];
end
%%
% figure
% plot(time,rad2deg(omega-omega(1,:)))
% title('Perigee Argument')
% ylabel('deg')
% legend(satLegend{:})
% grid on
% figure
% plot(time,rad2deg(i0-i0(1,:)))
% title('Ascending node inclination')
% ylabel('deg')
% legend(satLegend{:})
% grid on
% figure
% plot(time,sqrta-sqrta(1,:))
% title('Square root of major semi-axis')
% legend(satLegend{:})
% grid on
% ylabel('$\sqrt{m}$','interpreter','latex')
% save rinex_data1 -v7.3 
