function [satPos]=satAllPosition(nav,satId,time)
M=length(satId);
satPos=zeros(3,M);
%% Constants
mu=3986004.418e8;% gravitational constant (m^3/s^2)
we=7292115.0e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
for m=1:M
    id=find(nav(1,:)==satId(m));        %select all ephemeris for the satellite
    [~,j]=min(abs(time-nav(18,id))); %seek the more recent orbit parameters
    j=id(j);                            %select ephemeris to compute orbit
    toe=nav(18,j);    %time of ephemeris
    M0=nav(3,j);      %mean anomaly
    sqrta=nav(4,j);   %semi major-axis square-root
    a=sqrta^2;        %semi-major axis
    deltan=nav(5,j);  %mean motion difference
    ecc=nav(6,j);     %satellite orbit eccentricity
    omega=nav(7,j);   %argument of perigee
    cuc=nav(8,j);     %cos lat arg correction
    cus=nav(9,j);     %sin lat arg correction
    crc=nav(10,j);    %cos orbital radius correction
    crs=nav(11,j);    %sin orbitatl radius correction
    i0=nav(12,j);     %inclination at reference epoch (toe)
    doti0=nav(13,j);  %rate of inclination angle
    cic=nav(14,j);    %cos inclination correction
    cis=nav(15,j);    %sin inclination correction
    Omega=nav(16,j); %Ascending node's right ascension
    dotOmega=nav(17,j); %Rate of node's right ascension
    %%
    t=time-toe;
    if t>302400
        t=t-604800;
    elseif t<-302400
        t=t+604800;
    end
    M=M0+(sqrt(mu)/(sqrta*a)+deltan)*t;
    E=M;
    while abs(E-M-ecc*sin(E))>1e-15
        E=M+ecc*sin(E);
    end
    %true anomaly
    v=atan(sqrt(1-ecc^2)*sin(E)/(cos(E)-ecc));
    %orbit latitude (true anomaly+arg perigee)
    u=omega+v+cuc*cos(2*(omega+v))+cus*sin(2*(omega+v));
    r=a*(1-ecc*cos(E))+crc*cos(2*(omega+v))+crs*sin(2*(omega+v));
    i=i0+doti0*t+cic*cos(2*(omega+v))+cis*sin(2*(omega+v));
    Omega=Omega+(dotOmega-we)*t-we*toe;
    C=rotZ(-Omega)*rotX(-i)*rotZ(-u);
    satPos(:,m)=C*[r;0;0];
end
end