function skyPlot(obs,nav,p0,elevMask)
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
sats=unique(nav(1,:));
Nsats=numel(sats);
if Nsats>10
    rng(123);
    satColor=rand(Nsats,3);
end
satColor(1,:)=[0 0 1];
satColor(2,:)=[0 1 0];
satColor(3,:)=[1 0 0];
satColor(4,:)=[0 0.4470 0.7410];
satColor(5,:)=[0.8500 0.3250 0.0980];
satColor(6,:)=[0.9290 0.6940 0.1250];
satColor(7,:)=[0.4940 0.1840 0.5560];
satColor(8,:)=[0.4660 0.6740 0.1880];
satColor(9,:)=[0.3010 0.7450 0.9330];
satColor(10,:)=[0.6350 0.0780 0.1840];
time=unique(obs(:,2));
lla=SingleLlaFromEcef(p0);
figure
ax=polaraxes;
ax.RTick=0:15:90;
ax.RTickLabel=90:-15:0;
ax.RLim=[0 90];
hold on
N=length(time);
azim=zeros(Nsats,N);
elev=zeros(Nsats,N);
fprintf("Sky plot...  ")
satsOn=unique(obs(:,4));
for k=1:N
    epoch=obs(obs(:,2)==time(k),:);
    M=size(epoch,1);
    for m=1:M
        satId=epoch(m,4);
        id=find(nav(1,:)==satId);  %select all ephemeris for the satellite
        if size(id,2)>0
            [~,j]=min(abs(time(k)-nav(18,id)));  %seek for the most recent orbit parameters
            j=id(j); %nav msg for the m-th sat
            tems=time(k)-epoch(m,5)/c; %(GPST)
            psat=satPosition(nav(:,j),tems);
            theta=wie*norm(psat-p0)/c;
            psat=rotZ(theta)*psat;
            Cen=DCM_en(lla(1),lla(2));
            los=Cen'*(psat-p0)/norm(psat-p0); %Line of Sight (LOS)
            azim(sats==satId,k)=atan2d(los(1),los(2)); %Azimute (NED frame)
            elev(sats==satId,k)=asind(-los(3)); %Elevation (rad) (NED frame)
            %         polarplot(azim,90-elev,'o','color',satColor(sats==satId,:))
            %         drawnow
        end
    end
    if ~mod(k,round(N/10))
        fprintf('*');
    end
end
fprintf("\n")
%%
azim(azim==0)=NaN;
elev(elev==0)=NaN;
Non=numel(satsOn);
for m=1:Non
    p{m}=polarplot(deg2rad(azim(sats==satsOn(m),:)),90-elev(sats==satsOn(m),:),'.','color',satColor(m,:));
    strLeg{m}=['G' num2str(satsOn(m))];
    text(deg2rad(azim(sats==satsOn(m),end)), 90-elev(sats==satsOn(m),end), ['G' num2str(satsOn(m))], 'horiz', 'center', 'vert', 'top', 'rotation', 0)
end
theta=0:0.1:2*pi;
polarplot(theta,(90-elevMask)*ones(length(theta)),'k--')
legend([p{:}],strLeg)
title('Sky Plot')
end

