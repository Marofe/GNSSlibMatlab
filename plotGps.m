load('rinex_data1.mat')
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
filename='gpsOrbit.gif';
frame = getframe;
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
for i=1:360
    view(i,10)
        frame = getframe;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    %      Write to the GIF File
    imwrite(imind,cm,filename,'gif','WriteMode','append');
    pause(0.1)
end