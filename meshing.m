function [connec,nodecoord] = meshing(xdim,ydim,Ny,Nx)

%input:
%1. xdim: total x length (m)
%2. ydim: total y length (m)
%3. Ny: number of rectangular elements in y direction
%4. Nx: number of rectangular elements in x direction

%output:
%1. connec: connectivity matrix -> %[ele #, node1, node2, node3, node4]
%2. nodecoord: nodal coordinates matrix -> %[node #, x coord, y coord, z coord]

xlen = xdim/Nx;
ylen = ydim/Ny;

nodecoord = zeros((Ny+1)*(Nx+1),4); %[node #, x coord, y coord, z coord]
connec = zeros(Ny*Nx,5); %[ele #, node1, node2, node3, node4]

for i = 1:Nx+1 %loop over the number of "columns" (x division lines)
    nodecoord((i-1)*(Ny+1)+1:i*(Ny+1),1) = linspace((i-1)*(Ny+1)+1,i*(Ny+1),Ny+1); %coord number 
    nodecoord((i-1)*(Ny+1)+1:i*(Ny+1),2) = (i-1)*xlen; %x pos 
    nodecoord((i-1)*(Ny+1)+1:i*(Ny+1),3) = linspace(0,Ny,Ny+1)*ylen; %y pos
    nodecoord((i-1)*(Ny+1)+1:i*(Ny+1),4) = zeros(Ny+1,1); %z pos
    
    if i > 1 %assign 4 nodes to each rectangular element 
        connec((i-2)*Ny+1:(i-1)*Ny,1) = linspace((i-2)*Ny+1,(i-1)*Ny,Ny);
        connec((i-2)*Ny+1:(i-1)*Ny,2) = nodecoord((i-2)*(Ny+1)+1:1:(i-1)*(Ny+1)-1,1);
        connec((i-2)*Ny+1:(i-1)*Ny,3) = nodecoord((i-1)*(Ny+1)+1:1:(i)*(Ny+1)-1,1);
        connec((i-2)*Ny+1:(i-1)*Ny,4) = nodecoord((i-1)*(Ny+1)+2:1:(i)*(Ny+1),1);
        connec((i-2)*Ny+1:(i-1)*Ny,5) = nodecoord((i-2)*(Ny+1)+2:1:(i-1)*(Ny+1),1);
    end
end
% % 
figure,
%plotting the mesh
for i = 1:size(connec,1)
    plot([nodecoord(connec(i,2),2),nodecoord(connec(i,3),2)],[nodecoord(connec(i,2),3),nodecoord(connec(i,3),3)],'b');
    hold on,
    plot([nodecoord(connec(i,3),2),nodecoord(connec(i,4),2)],[nodecoord(connec(i,3),3),nodecoord(connec(i,4),3)],'b');
    hold on,
    plot([nodecoord(connec(i,4),2),nodecoord(connec(i,5),2)],[nodecoord(connec(i,4),3),nodecoord(connec(i,5),3)],'b');
    hold on,
    plot([nodecoord(connec(i,2),2),nodecoord(connec(i,5),2)],[nodecoord(connec(i,2),3),nodecoord(connec(i,5),3)],'b');
    hold on,
    text(0.5*(nodecoord(connec(i,2),2)+nodecoord(connec(i,3),2)),0.5*(nodecoord(connec(i,2),3)+nodecoord(connec(i,4),3)),int2str(connec(i,1)),'Color','red','FontSize',10);
    hold on,
    text(nodecoord(connec(i,2),2),nodecoord(connec(i,2),3),int2str(nodecoord(connec(i,2),1)),'FontSize',7);
    hold on,
    text(nodecoord(connec(i,3),2),nodecoord(connec(i,3),3),int2str(nodecoord(connec(i,3),1)),'FontSize',7);
    hold on,
    text(nodecoord(connec(i,4),2),nodecoord(connec(i,4),3),int2str(nodecoord(connec(i,4),1)),'FontSize',7);
    hold on,
    text(nodecoord(connec(i,5),2),nodecoord(connec(i,5),3),int2str(nodecoord(connec(i,5),1)),'FontSize',7);
    hold on,
end
end 
