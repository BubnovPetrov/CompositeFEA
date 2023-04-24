%% case non-conform CLT simple tensile 
clear all
clc 

xdim = 5;
ydim = 5;
Ny = 3; 
Nx = 3; 

BC.pos = {{0,0};
    %{0,ydim};
    {0,[0,ydim]};
    {[0,xdim],ydim};
    {xdim,0}};

BC.type = {{'x','y','w','wx','wy'};
    {'x'};
    %{'Nx'};
    {'w'};
    {'Nx'}};

BC.val = {{0,0,0,0,0};
    {0};
    {-100};
    {100}};

%% Square Plate Bending - RQNC4 Element
clear all
clc 

xdim = 3;
ydim = 3;

%mesh size 2x2: Ny = 2; Nx = 2; 
%mesh size 6x6: Ny = 4; Nx = 4; 
%mesh size 20x20: Ny = 20; Nx = 20;
%mesh size 30x30: 
Ny = 30; Nx = 30;
%mesh size 100x100: Ny = 100; Nx = 100;

BC.pos = {{[0,xdim],ydim};
    {[0,xdim],0};
    {xdim,[0,ydim]};
    {0,[0,ydim]};
    {[0,xdim],[0,ydim]}};

BC.type = {{'w','wx','wy'};
    {'w','wx','wy'};
    {'w','wx','wy'};
    {'w','wx','wy'};
    {'q'}};

BC.val = {{0,0,0};
    {0,0,0};
    {0,0,0};
    {0,0,0};
    {-5}};

%ply: [G12,E1,E2,v12,v21,theta,h]
%case1a): AR = 1; L/T = 100
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.01;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.01;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.01];

%case1b): AR = 1; L/T = 50
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.02;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.02;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.02];

%case1c): AR = 1; L/T = 25
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.04;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.04;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.04];

%case1d): AR = 1; L/T = 5
ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.2;
    3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.2;
    3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.2];

%case1d): AR = 1; L/T = 10
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.1;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.1;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.1];

[connec,nodecoord] = meshing(xdim,ydim,Ny,Nx);


%% Cantilever Test Case - RQNC4 Element
clear all
clc 
xdim = 10; ydim = 2;
%mesh size 10x2: Ny = 2; Nx = 10; 
%mesh size 20x4: Ny = 4; Nx = 20;
%mesh size 40x8: 
Ny = 8; Nx = 40;
%mesh size 80x16: Ny = 16; Nx = 80;
%mesh size 160x32: Ny = 32; Nx = 160;

BC.pos = {{xdim,[0,ydim]};
    {0,[0,ydim]}};

BC.type = {{'Q'};
    {'w','wx','x','y','wy'}};

BC.val = {{10};
    {0,0,0,0,0}};

E2 = 10^8;
%AR = 100 - [0/90/90/0]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005];

% %AR = 10 - [0/90/90/0]
ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.05;
    0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.05;
    0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.05;
    0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.05];

%AR = 100 - [0/90/0/90]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005];

% %AR = 10 - [0/90/0/90]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.05];

% %AR = 10 - [45/-45/45/-45]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.05];


% %AR = 100 - [45/-45/45/-45]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.005];

[connec,nodecoord] = meshing(xdim,ydim,Ny,Nx);

%% Nonconforming RBM tests 
% RIGID BODY MOTION FOR SINGLE ELEMENT NO LAYER

clear all
clc 
xdim = 1; ydim = 1;
%Nx = 1; Ny = 1; %1x1 mesh
Nx = 2; Ny = 2; %2x2 mesh
ply = [4.8*10^9,230*10^9,230*10^9,0.25,0.25,0,0.005];
[connec,nodecoord] = meshing(xdim,ydim,Ny,Nx);

%% Nonconforming element case without Von Karman Strain
ktot = sparse(5*max(nodecoord(:,1)),5*max(nodecoord(:,1))); %5 dof per node
p = 2; 
for i = 1:size(connec,1) %loop over elements
    x_e = nodecoord(connec(i,2:5),2); 
    y_e = nodecoord(connec(i,2:5),3); 
    [KA,KB,KC,KD] = elelaminatestiff_noncon(p,x_e,y_e,ply);
    
    rowpos1 = [connec(i,2)*2-1,connec(i,2)*2,connec(i,3)*2-1,connec(i,3)*2,connec(i,4)*2-1,connec(i,4)*2,connec(i,5)*2-1,connec(i,5)*2];
    rowpos2 = [connec(i,2)*3-2,connec(i,2)*3-1,connec(i,2)*3,connec(i,3)*3-2,connec(i,3)*3-1,connec(i,3)*3,...
            connec(i,4)*3-2,connec(i,4)*3-1,connec(i,4)*3,connec(i,5)*3-2,connec(i,5)*3-1,connec(i,5)*3]+ 2*max(nodecoord(:,1));
    for j = 1:8 %loop over the first 8 columns, populate rows
        colpos1 = (connec(i,ceil(j/2)+1)*2-mod(j,2))*ones(1,8);       
        colpos2 = (connec(i,ceil(j/2)+1)*2-mod(j,2))*ones(1,12); 
        tol = max(eps(ktot(rowpos1,colpos1(1))),eps(KA(:,j)));
        ktot = ktot + sparse(rowpos1,colpos1,KA(:,j),5*max(nodecoord(:,1)),5*max(nodecoord(:,1)));
        ktot(rowpos1(abs(ktot(rowpos1,colpos1(1)))<=100*tol),colpos1(1)) = 0;
        tol = max(eps(ktot(rowpos2,colpos2(1))),eps(KC(:,j)));
        ktot = ktot + sparse(rowpos2,colpos2,KC(:,j),5*max(nodecoord(:,1)),5*max(nodecoord(:,1)));
        ktot(rowpos2(abs(ktot(rowpos2,colpos2(1)))<=100*tol),colpos2(1)) = 0;
    end
    for j = 1:12 % loop over the last 12 columns
        colpos1 = (2*max(nodecoord(:,1))+connec(i,ceil(j/3)+1)*3-((ceil(j/3)-floor(j/3))*3-mod(j,3)))*ones(1,8);
        colpos2 = (2*max(nodecoord(:,1))+connec(i,ceil(j/3)+1)*3-((ceil(j/3)-floor(j/3))*3-mod(j,3)))*ones(1,12);
        tol = max(eps(ktot(rowpos1,colpos1(1))),eps(KB(:,j)));
        ktot = ktot + sparse(rowpos1,colpos1,KB(:,j),5*max(nodecoord(:,1)),5*max(nodecoord(:,1)));
        ktot(rowpos1(abs(ktot(rowpos1,colpos1(1)))<=100*tol),colpos1(1)) = 0;
        tol = max(eps(ktot(rowpos2,colpos2(1))),eps(KD(:,j)));
        ktot = ktot + sparse(rowpos2,colpos2,KD(:,j),5*max(nodecoord(:,1)),5*max(nodecoord(:,1)));
        ktot(rowpos2(abs(ktot(rowpos2,colpos2(1)))<=100*tol),colpos2(1)) = 0;
    end
end 

[d,d_dof,f] = solve_noncon(connec,nodecoord,BC,ktot,p);
ele = 1;
graphing(connec,nodecoord,d,ele,f)

%% Rigid body mode test case for nonconforming element - 1x1 mesh 
% Note the number of integration point for the Gauss quadrature needs to be
% 3 for the eigenvalues to not be imaginary.

K = full(ktot);
xmov = [repmat([1,0],1,4),zeros(1,20-8)];
ymov = [repmat([0,1],1,4),zeros(1,20-8)];
zmov = [zeros(1,8),repmat([1,0,0],1,4)];
f = zeros(length(zmov),1);

theta = 1;
zrot = [0,0,-15*sind(theta),-15+15*cosd(theta),...
    -(xdim-xdim*cosd(theta)),xdim*sind(theta),...
    -(xdim-(xdim*cosd(theta)-xdim*sind(theta))), xdim*sind(theta)+xdim*cosd(theta)-xdim,...
    zeros(1,12)];

[q,v] = eigs(K,6,"smallestab");

kk = q(1:8,4:6);
t1 = repmat([0,1],1,4);
t2 = repmat([1,0],1,4);
s1 = kk\t1';
s2 = kk\t2';
s3 = null([s1,s2]');
c3 = kk*s3;
graphing(connec,nodecoord,[c3./10^7.5;zeros(12,1)],1,f);

kz = q(9:end,1:3);
w1 = repmat([1,0,0],1,4);
ss1 = kz\w1';
sz = null([ss1]');
cz = kz*sz;
rx = cz*null(cz(2,:));
ry = cz*null(cz(3,:));

graphing(connec,nodecoord,[zeros(8,1);ry./10^2],1,f);
graphing(connec,nodecoord,[zeros(8,1);rx./10^2],1,f);

%% Rigid body mode test case for nonconforming element - 2x2 mesh 
% Note the number of integration point for the Gauss quadrature needs to be
% 3 for the eigenvalues to not be imaginary.

K = full(ktot);
xmov = [repmat([1,0],1,9),zeros(1,45-18)];
ymov = [repmat([0,1],1,9),zeros(1,45-18)];
zmov = [zeros(1,18),repmat([1,0,0],1,9)];
f = zeros(length(zmov),1);
theta = 1;
zrot = [0,0,-15*sind(theta),-15+15*cosd(theta),...
    -(xdim-xdim*cosd(theta)),xdim*sind(theta),...
    -(xdim-(xdim*cosd(theta)-xdim*sind(theta))), xdim*sind(theta)+xdim*cosd(theta)-xdim,...
    zeros(1,12)];

[q,v] = eigs(K,6,"smallestab");

kk = q(1:18,4:6);
t1 = repmat([0,1],1,9);
t2 = repmat([1,0],1,9);
s1 = kk\t1';
s2 = kk\t2';
s3 = null([s1,s2]');
c3 = kk*s3;
graphing(connec,nodecoord,[c3./10^7.5;zeros(45-18,1)],1,f);

kz = q(19:end,1:3);
w1 = repmat([1,0,0],1,9);
ss1 = kz\w1';
sz = null([ss1]');
cz = kz*sz;
rx = cz*null(cz(2,:));
ry = cz*null(cz(3,:));

graphing(connec,nodecoord,[zeros(18,1);ry./10^2],1,f);
graphing(connec,nodecoord,[zeros(18,1);rx./10^2],1,f);

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Square Plate Case - RQC4 Element
clear all
clc 
xdim = 3; ydim = 3;
%mesh size 2x2: Ny = 2; Nx = 2; 
%mesh size 6x6: Ny = 6; Nx = 6; 
%mesh size 20x20: Ny = 20; Nx = 20;
%mesh size 30x30: 
Ny = 30; Nx = 30;
%mesh size 100x100: Ny = 100; Nx = 100;
%mesh size 150x150: Ny = 150; Nx = 150;

BC.pos = {{[0,xdim],ydim};
    {[0,xdim],0};
    {xdim,[0,ydim]};
    {0,[0,ydim]};
    {[0,xdim],[0,ydim]}};

BC.type = {{'w','wx','wy','wxy'};
    {'w','wx','wy','wxy'};
    {'w','wx','wy','wxy'};
    {'w','wx','wy','wxy'};
    {'q'}};

BC.val = {{0,0,0,0};
    {0,0,0,0};
    {0,0,0,0};
    {0,0,0,0};
    {-5}};
    
%ply: [G12,E1,E2,v12,v21,theta,h]
%case1a): AR = 1; L/T = 100
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.01;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.01;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.01];

%case1b): AR = 1; L/T = 50
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.02;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.02;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.02];

%case1c): AR = 1; L/T = 25
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.04;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.04;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.04];

%case1d): AR = 1; L/T = 5
ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.2;
    3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.2;
    3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.2];

%case1d): AR = 1; L/T = 10
% ply = [3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.1;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,pi/2,0.1;
%     3447.38*10^6,172369*10^6,6894.76*10^6,0.25,6894.76*10^6*0.25/172369*10^6,0,0.1];

[connec,nodecoord] = meshing(xdim,ydim,Ny,Nx);

%% Cantilever Test Case - RQC4 Element
clear all
clc 
xdim = 10; ydim = 2;

%mesh size 10x2: 
Ny = 2; Nx = 10; 
%mesh size 20x4: Ny = 4; Nx = 20;
%mesh size 40x8: Ny = 8; Nx = 40;
%mesh size 80x16: Ny = 16; Nx = 80;
%mesh size 160x32: Ny = 32; Nx = 160;

BC.pos = {{xdim,[0,ydim]};
    {0,[0,ydim]}};

BC.type = {{'Q'};
    {'w','wx','x','y','wy','wxy'}};

BC.val = {{10};
    {0,0,0,0,0,0}};

E2 = 10^8;
%AR = 100 [0/90/90/0]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005];


% %AR = 10 [0/90/90/0]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(0.5*E2),0,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(0.5*E2),pi/2,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(0.5*E2),pi/2,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(0.5*E2),0,0.05];

%AR = 100 -[0/90/0/90]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.005;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.005];

% %AR = 10 - [0/90/0/90]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),0,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/2,0.05];

% %AR = 10 - [45/-45/45/-45]
% ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.05;
%     0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.05];

% %AR = 100 - [45/-45/45/-45]
ply = [0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.005;
    0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.005;
    0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),pi/4,0.005;
    0.5*E2,25*E2,E2,0.25,E2*0.25/(25*E2),-pi/4,0.005];

[connec,nodecoord] = meshing(xdim,ydim,Ny,Nx);

%%
% RIGID BODY MOTION FOR SINGLE ELEMENT NO LAYER

clear all
clc 
xdim = 1; ydim = 1;
Nx = 1; Ny = 1;
%Nx = 2; Ny = 2;
ply = [4.8*10^9,230*10^9,230*10^9,0.25,0.25,0,0.005];
%ply = [2*10^9,230*10^9,230*10^9,0.25,0.25,0,0.005];
[connec,nodecoord] = meshing(xdim,ydim,Ny,Nx);

%% Conforming element case without Von Karman Strain
ktot = sparse(6*max(nodecoord(:,1)),6*max(nodecoord(:,1))); %6 dof per node
p = 3; 
for i = 1:size(connec,1) %loop over elements
    x_e = nodecoord(connec(i,2:5),2); 
    y_e = nodecoord(connec(i,2:5),3); 
    [KA,KB,KC,KD] = elelaminatestiff_con(p,x_e,y_e,ply);
    
    rowpos1 = [connec(i,2)*2-1,connec(i,2)*2,connec(i,3)*2-1,connec(i,3)*2,connec(i,4)*2-1,connec(i,4)*2,connec(i,5)*2-1,connec(i,5)*2];
    rowpos2 = [connec(i,2)*4-3,connec(i,2)*4-2,connec(i,2)*4-1,connec(i,2)*4,connec(i,3)*4-3,connec(i,3)*4-2,connec(i,3)*4-1,connec(i,3)*4,...
            connec(i,4)*4-3,connec(i,4)*4-2,connec(i,4)*4-1,connec(i,4)*4,connec(i,5)*4-3,connec(i,5)*4-2,connec(i,5)*4-1,connec(i,5)*4]+ 2*max(nodecoord(:,1));
    for j = 1:8 %loop over the first 8 columns, populate rows
        colpos1 = (connec(i,ceil(j/2)+1)*2-mod(j,2))*ones(1,8);       
        colpos2 = (connec(i,ceil(j/2)+1)*2-mod(j,2))*ones(1,16); 
   
        tol = max(eps(ktot(rowpos1,colpos1(1))),eps(KA(:,j)));
        ktot = ktot + sparse(rowpos1,colpos1,KA(:,j),6*max(nodecoord(:,1)),6*max(nodecoord(:,1)));
        ktot(rowpos1(abs(ktot(rowpos1,colpos1(1)))<=100*tol),colpos1(1)) = 0;
        tol = max(eps(ktot(rowpos2,colpos2(1))),eps(KC(:,j)));
        ktot = ktot + sparse(rowpos2,colpos2,KC(:,j),6*max(nodecoord(:,1)),6*max(nodecoord(:,1)));
        ktot(rowpos2(abs(ktot(rowpos2,colpos2(1)))<=100*tol),colpos2(1)) = 0;
    end
    for j = 1:16 % loop over the last 16 columns
        colpos1 = (2*max(nodecoord(:,1))+connec(i,ceil(j/4)+1)*4-((ceil(j/4)-floor(j/4))*4-mod(j,4)))*ones(1,8);
        colpos2 = (2*max(nodecoord(:,1))+connec(i,ceil(j/4)+1)*4-((ceil(j/4)-floor(j/4))*4-mod(j,4)))*ones(1,16);
        tol = max(eps(ktot(rowpos1,colpos1(1))),eps(KB(:,j)));
        ktot = ktot + sparse(rowpos1,colpos1,KB(:,j),6*max(nodecoord(:,1)),6*max(nodecoord(:,1)));
        ktot(rowpos1(abs(ktot(rowpos1,colpos1(1)))<=100*tol),colpos1(1)) = 0;
        tol = max(eps(ktot(rowpos2,colpos2(1))),eps(KD(:,j)));
        ktot = ktot + sparse(rowpos2,colpos2,KD(:,j),6*max(nodecoord(:,1)),6*max(nodecoord(:,1)));
        ktot(rowpos2(abs(ktot(rowpos2,colpos2(1)))<=100*tol),colpos2(1)) = 0;
    end
end 

[d,d_dof,f] = solve_con(connec,nodecoord,BC,ktot,p);
ele = 2;
graphing(connec,nodecoord,d,ele,f)

%% Rigid body mode test case for conforming element - 1x1 mesh
K = full(ktot);
%1x1
xmov = [repmat([1,0],1,4),zeros(1,24-8)];
ymov = [repmat([0,1],1,4),zeros(1,24-8)];
zmov = [zeros(1,8),repmat([1,0,0,0],1,4)];
f = zeros(length(zmov),1);
[q,v] = eigs(K,7,'smallestab');
%qt = [q(:,3:5),q(:,1:2),q(:,6)];
qt = [q(:,1:3),q(:,4:6)];

kk = qt(1:8,4:6);
t1 = [1,0,1,0,1,0,1,0];
t2 = [0,1,0,1,0,1,0,1];
s1 = kk\t1';
s2 = kk\t2';
s3 = null([s1,s2]');
c3 = kk*s3;
graphing(connec,nodecoord,[c3./10^7.5;zeros(16,1)],2,f);

pos1 = [0,0]' + c3(1:2);
pos2 = [0,ydim]' + c3(3:4);
pos3 = [xdim,0]' + c3(5:6);
pos4 = [xdim,ydim]' + c3(7:8);
l1 = norm(pos2-pos1);
l2 = norm(pos3-pos1);
ll1 = norm(pos4-pos3);
ll2 = norm(pos4-pos2);

kz = qt(9:24,1:3);
w1 = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];
ss1 = kz\w1';
sz = null([ss1]');
cz = kz*sz;
rx = cz*null(cz(2,:));
ry = cz*null(cz(3,:));


%1x1 mesh, p=2, "hourglassing/zero energy RBM mode":
% [q,v] = eigs(K,7,'smallestab');
% qt = [q(:,1:3),q(:,4:7)];
% kz = qt(9:24,1:4);
% w1 = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];
% ss1 = kz\w1';
% sz = null([ss1]');
% cz = kz*sz;
% rx = cz*null(cz(1:2,:));
% ry = cz*null(cz(2:3,:));
% rxy_ss1 = null([cz\rx,cz\ry]');
% rxy = cz*rxy_ss1;
%after this, run the troubleshooting_wxy.m code to get the "RBM shape" and
%the corresponding strain values at guass points which shows that all the
%strain values are near zero, thus making this a zero energy mode

graphing(connec,nodecoord,[zeros(8,1);rx./10^2],2,f);
graphing(connec,nodecoord,[zeros(8,1);ry./10^2],2,f);

theta = 5;
yrot = [zeros(1,4),-xdim*(1-cosd(theta)),0,-xdim*(1-cosd(theta)),0,zeros(1,6),1,0,0,1,0,0];

%% Rigid body mode test case for conforming element - 2x2 mesh
K = full(ktot);
xmov = [repmat([1,0],1,9),zeros(1,54-18)];
ymov = [repmat([0,1],1,9),zeros(1,54-18)];
zmov = [zeros(1,18),repmat([1,0,0,0],1,9)];
f = zeros(length(zmov),1);

[q,v] = eigs(K,6,'smallestab');
%qt = [q(:,3:5),q(:,1:2),q(:,6)];
qt = [q(:,1:3),q(:,4:6)];

kk = qt(1:18,4:6);
t1 = [repmat([1,0],1,9)];
t2 = [repmat([0,1],1,9)];
s1 = kk\t1';
s2 = kk\t2';
s3 = null([s1,s2]');
c3 = kk*s3;
graphing(connec,nodecoord,[c3;zeros(54-18,1)],2,f);

kz = qt(19:end,1:3);
w1 = [repmat([1,0,0,0],1,9)];
ss1 = kz\w1';
sz = null([ss1]');
cz = kz*sz;
rx = cz*null(cz(2,:));
ry = cz*null(cz(3,:));
%rxy_ss1 = null([cz\rx,cz\ry]');
%rxy = cz*rxy_ss1;
graphing(connec,nodecoord,[zeros(18,1);rx],2,f);
graphing(connec,nodecoord,[zeros(18,1);ry],2,f);


