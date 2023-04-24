function graphing(connec,nodecoord,d,ele,f)
%Input variable:
%ele: 1- non-conforming element
%     2- conforming element

% plotting code deformed surface contour 
fh(1) = figure(1);
title('Deformed Mesh')
ax(1) = axes('parent',fh(1));
fh(2) = figure(2);
title('U Displacement Map')
ax(2) = axes('parent',fh(2));
fh(3) = figure(3);
title('V Displacement Map')
ax(3) = axes('parent',fh(3));
fh(4) = figure(4);
title('W Displacement Map')
ax(4) = axes('parent',fh(4));
figure,
for i = 1:size(connec,1) %loop over elements
    x_e = nodecoord(connec(i,2:5),2); 
    y_e = nodecoord(connec(i,2:5),3); 
    
    xtrans = @(x,x_e) (2*x-x_e(1)-x_e(2))/(x_e(2)-x_e(1)); %this is a
    ytrans = @(y,y_e) (2*y-y_e(1)-y_e(4))/(y_e(4)-y_e(1)); %this is b 
    
    %Bxy is a 2 by 8 matrix (in-plane shape functions)
    Bxy = @(a,b)reshape([(a./4.0-1.0./4.0).*(b-1.0),0.0,0.0,(a./4.0-1.0./4.0).*(b-1.0),-(a./4.0+1.0./4.0).*(b-1.0),0.0,0.0, ...
        -(a./4.0+1.0./4.0).*(b-1.0),(a./4.0+1.0./4.0).*(b+1.0),0.0,0.0,(a./4.0+1.0./4.0).*(b+1.0),-(a./4.0-1.0./4.0).*(b+1.0), ...
        0.0,0.0,-(a./4.0-1.0./4.0).*(b+1.0)],[2,8]);

    if ele == 1
        %Bw is a 1 by 12 matrix
        Bw = @(a,b)[-(a./8.0-1.0./8.0).*(b-1.0).*(a+b+a.^2+b.^2-2.0),(a-1.0).^2.*(a+1.0).*(b-1.0).*(-1.0./8.0), ...
            (a-1.0).*(b-1.0).^2.*(b+1.0).*(-1.0./8.0),(a./8.0+1.0./8.0).*(b-1.0).*(-a+b+a.^2+b.^2-2.0), ...
            (a-1.0).*(a+1.0).^2.*(b-1.0).*(-1.0./8.0),((a+1.0).*(b-1.0).^2.*(b+1.0))./8.0, ...
            (a./8.0+1.0./8.0).*(b+1.0).*(a+b-a.^2-b.^2+2.0),((a-1.0).*(a+1.0).^2.*(b+1.0))./8.0, ...
            ((a+1.0).*(b-1.0).*(b+1.0).^2)./8.0,(a./8.0-1.0./8.0).*(b+1.0).*(a-b+a.^2+b.^2-2.0), ...
            ((a-1.0).^2.*(a+1.0).*(b+1.0))./8.0,(a-1.0).*(b-1.0).*(b+1.0).^2.*(-1.0./8.0)];
        
    elseif ele == 2
    %Bw is a 1 by 16 matrix
        Bw = @(a,b)[((a-1.0).^2.*(a+2.0).*(b-1.0).^2.*(b+2.0))./1.6e+1,((a-1.0).^2.*(a+1.0).*(b-1.0).^2.*(b+2.0))./1.6e+1, ...
            ((a-1.0).^2.*(a+2.0).*(b-1.0).^2.*(b+1.0))./1.6e+1,((a-1.0).^2.*(a+1.0).*(b-1.0).^2.*(b+1.0))./1.6e+1, ...
            (a+1.0).^2.*(a-2.0).*(b-1.0).^2.*(b+2.0).*(-1.0./1.6e+1),((a-1.0).*(a+1.0).^2.*(b-1.0).^2.*(b+2.0))./1.6e+1, ...
            (a+1.0).^2.*(a-2.0).*(b-1.0).^2.*(b+1.0).*(-1.0./1.6e+1),((a-1.0).*(a+1.0).^2.*(b-1.0).^2.*(b+1.0))./1.6e+1, ...
            ((a+1.0).^2.*(a-2.0).*(b+1.0).^2.*(b-2.0))./1.6e+1,(a-1.0).*(a+1.0).^2.*(b+1.0).^2.*(b-2.0).*(-1.0./1.6e+1), ...
            (a+1.0).^2.*(a-2.0).*(b-1.0).*(b+1.0).^2.*(-1.0./1.6e+1),((a-1.0).*(a+1.0).^2.*(b-1.0).*(b+1.0).^2)./1.6e+1, ...
            (a-1.0).^2.*(a+2.0).*(b+1.0).^2.*(b-2.0).*(-1.0./1.6e+1),(a-1.0).^2.*(a+1.0).*(b+1.0).^2.*(b-2.0).*(-1.0./1.6e+1), ...
            ((a-1.0).^2.*(a+2.0).*(b-1.0).*(b+1.0).^2)./1.6e+1,((a-1.0).^2.*(a+1.0).*(b-1.0).*(b+1.0).^2)./1.6e+1];
    end
        
    xbin_s = linspace(x_e(1),x_e(2),100);
    ybin_s = linspace(y_e(1),y_e(4),100);
    [xbin,ybin] = meshgrid(xbin_s,ybin_s);
    
    curelem = connec(i,2:5);
    for m = 1:length(ybin_s)
        for n = 1:length(xbin_s)
            xydef = Bxy(xtrans(xbin_s(n),x_e),ytrans(ybin_s(m),y_e));
            xydef = xydef*[d(curelem(1)*2-1:curelem(1)*2);d(curelem(2)*2-1:curelem(2)*2);...
                        d(curelem(3)*2-1:curelem(3)*2);d(curelem(4)*2-1:curelem(4)*2)];
            xdefplot(n,m) = xydef(1);
            ydefplot(n,m) = xydef(2);
            
            zdef = Bw(xtrans(xbin_s(n),x_e),ytrans(ybin_s(m),y_e));
            if ele == 1
                zdefplot(n,m) = zdef*[d(2*max(nodecoord(:,1))+curelem(1)*3-2:2*max(nodecoord(:,1))+curelem(1)*3);...
                    d(2*max(nodecoord(:,1))+curelem(2)*3-2:2*max(nodecoord(:,1))+curelem(2)*3);...
                    d(2*max(nodecoord(:,1))+curelem(3)*3-2:2*max(nodecoord(:,1))+curelem(3)*3);...
                    d(2*max(nodecoord(:,1))+curelem(4)*3-2:2*max(nodecoord(:,1))+curelem(4)*3)];
            elseif ele == 2
                zdefplot(n,m) = zdef*[d(2*max(nodecoord(:,1))+curelem(1)*4-3:2*max(nodecoord(:,1))+curelem(1)*4);...
                    d(2*max(nodecoord(:,1))+curelem(2)*4-3:2*max(nodecoord(:,1))+curelem(2)*4);...
                    d(2*max(nodecoord(:,1))+curelem(3)*4-3:2*max(nodecoord(:,1))+curelem(3)*4);...
                    d(2*max(nodecoord(:,1))+curelem(4)*4-3:2*max(nodecoord(:,1))+curelem(4)*4)];
            end
            magdefplot(n,m) = sqrt(xydef(1)^2 + xydef(2)^2 + zdefplot(n,m)^2);
        end
    end
    xdefplotd = xdefplot + xbin;
    ydefplotd = ydefplot + ybin;
    xdefplot = subs(xydef(1,:),{x,y},{xbin,ybin})+xbin;
    ydefplot = subs(xydef(2,:),{x,y},{xbin,ybin})+ybin;
    zdefplot = subs(zdef,{x,y},{xbin,ybin});

    surf(double(xdefplotd)',double(ydefplotd)',double(zdefplot),'EdgeColor','none','Parent',ax(1)),
    hold on,
    %hold(ax(1),'on'),
    surf(xbin,ybin,xdefplot','EdgeColor','none','Parent',ax(2)),
    hold on,
    %hold(ax(2),'on'),
    surf(xbin,ybin,ydefplot','EdgeColor','none','Parent',ax(3)),
    hold on,
    %hold(ax(3),'on'),
    surf(xbin,ybin,zdefplot','EdgeColor','none','Parent',ax(4)),
    hold on,
    %hold(ax(4),'on'),
    
    %use the part below if you only want to plot the magnitude of
    %deformation. You must comment all the above graphing codes.
    
    %surf(xbin,ybin,magdefplot','EdgeColor','none'),
    %hold on,
end
set(fh(1),'Name','Deformed Mesh');
set(ax(2), 'xTickLabel',' ')
set(ax(2), 'yTickLabel',' ')
set(ax(2), 'zTickLabel',' ')

set(ax(3), 'xTickLabel',' ')
set(ax(3), 'yTickLabel',' ')
set(ax(3), 'zTickLabel',' ')

set(ax(4), 'xTickLabel',' ')
set(ax(4), 'yTickLabel',' ')
set(ax(4), 'zTickLabel',' ')

%use the part below if you only want to label the plot of the magnitude of
%deformation. You must comment all the above graph labelling codes.
%title('Displacement Magnitude Map')
%title('z Displacement Map')
%xlabel('X-Position')
%ylabel('Y-Position')


% plotting code reaction force contour 
fh(5) = figure(5);
title('X Reaction Force')
ax(5) = axes('parent',fh(5));
fh(6) = figure(6);
title('Y Reaction Force')
ax(6) = axes('parent',fh(6));
fh(7) = figure(7);
title('Z Reaction Force')
ax(7) = axes('parent',fh(7));

for i = 1:size(connec,1) %loop over elements
    x_e = nodecoord(connec(i,2:5),2); 
    y_e = nodecoord(connec(i,2:5),3); 
    
    xtrans = @(x,x_e) (2*x-x_e(1)-x_e(2))/(x_e(2)-x_e(1)); %this is a
    ytrans = @(y,y_e) (2*y-y_e(1)-y_e(4))/(y_e(4)-y_e(1)); %this is b 
    
    %Bxy is a 2 by 8 matrix (in-plane shape functions)
    Bxy = @(a,b)reshape([(a./4.0-1.0./4.0).*(b-1.0),0.0,0.0,(a./4.0-1.0./4.0).*(b-1.0),-(a./4.0+1.0./4.0).*(b-1.0),0.0,0.0, ...
        -(a./4.0+1.0./4.0).*(b-1.0),(a./4.0+1.0./4.0).*(b+1.0),0.0,0.0,(a./4.0+1.0./4.0).*(b+1.0),-(a./4.0-1.0./4.0).*(b+1.0), ...
        0.0,0.0,-(a./4.0-1.0./4.0).*(b+1.0)],[2,8]);

    if ele == 1
        %Bw is a 1 by 12 matrix
        Bw = @(a,b)[-(a./8.0-1.0./8.0).*(b-1.0).*(a+b+a.^2+b.^2-2.0),(a-1.0).^2.*(a+1.0).*(b-1.0).*(-1.0./8.0), ...
            (a-1.0).*(b-1.0).^2.*(b+1.0).*(-1.0./8.0),(a./8.0+1.0./8.0).*(b-1.0).*(-a+b+a.^2+b.^2-2.0), ...
            (a-1.0).*(a+1.0).^2.*(b-1.0).*(-1.0./8.0),((a+1.0).*(b-1.0).^2.*(b+1.0))./8.0, ...
            (a./8.0+1.0./8.0).*(b+1.0).*(a+b-a.^2-b.^2+2.0),((a-1.0).*(a+1.0).^2.*(b+1.0))./8.0, ...
            ((a+1.0).*(b-1.0).*(b+1.0).^2)./8.0,(a./8.0-1.0./8.0).*(b+1.0).*(a-b+a.^2+b.^2-2.0), ...
            ((a-1.0).^2.*(a+1.0).*(b+1.0))./8.0,(a-1.0).*(b-1.0).*(b+1.0).^2.*(-1.0./8.0)];
        
    elseif ele == 2
    %Bw is a 1 by 16 matrix
        Bw = @(a,b)[((a-1.0).^2.*(a+2.0).*(b-1.0).^2.*(b+2.0))./1.6e+1,((a-1.0).^2.*(a+1.0).*(b-1.0).^2.*(b+2.0))./1.6e+1, ...
            ((a-1.0).^2.*(a+2.0).*(b-1.0).^2.*(b+1.0))./1.6e+1,((a-1.0).^2.*(a+1.0).*(b-1.0).^2.*(b+1.0))./1.6e+1, ...
            (a+1.0).^2.*(a-2.0).*(b-1.0).^2.*(b+2.0).*(-1.0./1.6e+1),((a-1.0).*(a+1.0).^2.*(b-1.0).^2.*(b+2.0))./1.6e+1, ...
            (a+1.0).^2.*(a-2.0).*(b-1.0).^2.*(b+1.0).*(-1.0./1.6e+1),((a-1.0).*(a+1.0).^2.*(b-1.0).^2.*(b+1.0))./1.6e+1, ...
            ((a+1.0).^2.*(a-2.0).*(b+1.0).^2.*(b-2.0))./1.6e+1,(a-1.0).*(a+1.0).^2.*(b+1.0).^2.*(b-2.0).*(-1.0./1.6e+1), ...
            (a+1.0).^2.*(a-2.0).*(b-1.0).*(b+1.0).^2.*(-1.0./1.6e+1),((a-1.0).*(a+1.0).^2.*(b-1.0).*(b+1.0).^2)./1.6e+1, ...
            (a-1.0).^2.*(a+2.0).*(b+1.0).^2.*(b-2.0).*(-1.0./1.6e+1),(a-1.0).^2.*(a+1.0).*(b+1.0).^2.*(b-2.0).*(-1.0./1.6e+1), ...
            ((a-1.0).^2.*(a+2.0).*(b-1.0).*(b+1.0).^2)./1.6e+1,((a-1.0).^2.*(a+1.0).*(b-1.0).*(b+1.0).^2)./1.6e+1];
    end
        
    xbin_s = linspace(x_e(1),x_e(2),100);
    ybin_s = linspace(y_e(1),y_e(4),100);
    [xbin,ybin] = meshgrid(xbin_s,ybin_s);
    
    curelem = connec(i,2:5);
    for m = 1:length(ybin_s)
        for n = 1:length(xbin_s)
            xyfrea = Bxy(xtrans(xbin_s(n),x_e),ytrans(ybin_s(m),y_e));
            xyfrea = xyfrea*[f(curelem(1)*2-1:curelem(1)*2);f(curelem(2)*2-1:curelem(2)*2);...
                        f(curelem(3)*2-1:curelem(3)*2);f(curelem(4)*2-1:curelem(4)*2)];
            xfreaplot(n,m) = xyfrea(1);
            yfreaplot(n,m) = xyfrea(2);
            
            zfrea = Bw(xtrans(xbin_s(n),x_e),ytrans(ybin_s(m),y_e));
            if ele == 1
                zfreaplot(n,m) = zfrea*[f(2*max(nodecoord(:,1))+curelem(1)*3-2:2*max(nodecoord(:,1))+curelem(1)*3);...
                    f(2*max(nodecoord(:,1))+curelem(2)*3-2:2*max(nodecoord(:,1))+curelem(2)*3);...
                    f(2*max(nodecoord(:,1))+curelem(3)*3-2:2*max(nodecoord(:,1))+curelem(3)*3);...
                    f(2*max(nodecoord(:,1))+curelem(4)*3-2:2*max(nodecoord(:,1))+curelem(4)*3)].*10^2;
            elseif ele == 2
                zfreaplot(n,m) = zfrea*[f(2*max(nodecoord(:,1))+curelem(1)*4-3:2*max(nodecoord(:,1))+curelem(1)*4);...
                    f(2*max(nodecoord(:,1))+curelem(2)*4-3:2*max(nodecoord(:,1))+curelem(2)*4);...
                    f(2*max(nodecoord(:,1))+curelem(3)*4-3:2*max(nodecoord(:,1))+curelem(3)*4);...
                    f(2*max(nodecoord(:,1))+curelem(4)*4-3:2*max(nodecoord(:,1))+curelem(4)*4)].*10^2;
            end
            
        end
    end
    surf(xbin,ybin,xfreaplot','EdgeColor','none','Parent',ax(5)),
    hold on,
    hold(ax(5),'on'),
    surf(xbin,ybin,yfreaplot','EdgeColor','none','Parent',ax(6)),
    hold on,
    hold(ax(6),'on'),
    surf(xbin,ybin,zfreaplot','EdgeColor','none','Parent',ax(7)),
    hold on,
    hold(ax(7),'on'),
end
set(ax(5), 'xTickLabel',' ')
set(ax(5), 'yTickLabel',' ')
set(ax(5), 'zTickLabel',' ')

set(ax(6), 'xTickLabel',' ')
set(ax(6), 'yTickLabel',' ')
set(ax(6), 'zTickLabel',' ')

set(ax(7), 'xTickLabel',' ')
set(ax(7), 'yTickLabel',' ')
set(ax(7), 'zTickLabel',' ')

% plotting codes node-wise
newcoord = zeros(size(nodecoord)); 
newcoord(:,1) = nodecoord(:,1); 
for i = 1:size(nodecoord,1)
%change the "1" value to whatever scale you desire.

newcoord(i,2) = 1*d(2*i-1)+nodecoord(i,2); %xcoord 
newcoord(i,3) = 1*d(2*i)+nodecoord(i,3); %ycoord
if ele == 1
    newcoord(i,4) = 1*d(2*max(nodecoord(:,1))+3*i-2); %zcoord
elseif ele == 2
    newcoord(i,4) = 1*d(2*max(nodecoord(:,1))+4*i-3); %zcoord
end
end 
%plotting the deformed mesh and boundary force 
figure,
for i = 1:size(connec,1)
    plot3([newcoord(connec(i,2),2),newcoord(connec(i,3),2)],[newcoord(connec(i,2),3),newcoord(connec(i,3),3)],[newcoord(connec(i,2),4),newcoord(connec(i,3),4)],'b');
    hold on,
    plot3([newcoord(connec(i,3),2),newcoord(connec(i,4),2)],[newcoord(connec(i,3),3),newcoord(connec(i,4),3)],[newcoord(connec(i,3),4),newcoord(connec(i,4),4)],'b','HandleVisibility','off');
    hold on,
    plot3([newcoord(connec(i,4),2),newcoord(connec(i,5),2)],[newcoord(connec(i,4),3),newcoord(connec(i,5),3)],[newcoord(connec(i,4),4),newcoord(connec(i,5),4)],'b','HandleVisibility','off');
    hold on,
    plot3([newcoord(connec(i,2),2),newcoord(connec(i,5),2)],[newcoord(connec(i,2),3),newcoord(connec(i,5),3)],[newcoord(connec(i,2),4),newcoord(connec(i,5),4)],'b','HandleVisibility','off');
    
    hold on,
    plot([nodecoord(connec(i,2),2),nodecoord(connec(i,3),2)],[nodecoord(connec(i,2),3),nodecoord(connec(i,3),3)],'r');
    hold on,
    plot([nodecoord(connec(i,3),2),nodecoord(connec(i,4),2)],[nodecoord(connec(i,3),3),nodecoord(connec(i,4),3)],'r');
    hold on,
    plot([nodecoord(connec(i,4),2),nodecoord(connec(i,5),2)],[nodecoord(connec(i,4),3),nodecoord(connec(i,5),3)],'r');
    hold on,
    plot([nodecoord(connec(i,2),2),nodecoord(connec(i,5),2)],[nodecoord(connec(i,2),3),nodecoord(connec(i,5),3)],'r');
end
legend('Deformed Mesh','Undeformed Mesh','Location','Best')
title('Mesh Nodal Displacement Map')
xlabel('X-Position')
ylabel('Y-Position')
zlabel('Z-Position')
    
end 