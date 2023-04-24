function [d,d_dof,f] = solve_noncon(connec,nodecoord,BC,ktot,p)
%for non-conforming case
%BC is a MATLAB data struct that contains the boundary conditions applied,
%including the physical location, the type of BC, and the corresponding
%value or function handle.

%BC.type condition cases: x,y,w,wx,wy - Displacement BCs
    %                     Nx,Ny,Nxy - Distributed Axial Load BCs (N/m)
    %                     Fx,Fy,Fxy - Point Axial Load BCs (N)
    %                     Mx,My,Mxy - Distributed Axial Moment BCs (Nm/m)
    %                     q - Domain Transverse Load (N/m^2)
    %                     Q - Distributed Transverse Load BC(transverse shear force) (N/m)
    
d = zeros(5*max(nodecoord(:,1)),1);
f = zeros(5*max(nodecoord(:,1)),1);
d_dof = zeros(5*max(nodecoord(:,1)),1);

%loop through all BC input and construct the F (force) vector 
for i = 1:size(BC.pos,1)
    %identifying relevant nodes for displacement BCs
    if length(BC.pos{i}{1}) == 2 %range of coordinate specified for x-axis
        if length(BC.pos{i}{2}) == 2 %an area specified
            nodes = nodecoord(nodecoord(:,2)>=BC.pos{i}{1}(1)&nodecoord(:,2)<=BC.pos{i}{1}(2)& ...
                nodecoord(:,3)>=BC.pos{i}{2}(1)&nodecoord(:,3)<=BC.pos{i}{2}(2),1);
        else %only a line of x-coords
            nodes = nodecoord(nodecoord(:,2)>=BC.pos{i}{1}(1)&nodecoord(:,2)<=BC.pos{i}{1}(2)&BC.pos{i}{2}==nodecoord(:,3),1);
        end
    elseif length(BC.pos{i}{2}) == 2 %range of coordinate specified for y-axis (only a line of y-coord)
        nodes = nodecoord(nodecoord(:,3)>=BC.pos{i}{2}(1)&nodecoord(:,3)<=BC.pos{i}{2}(2)&BC.pos{i}{1}==nodecoord(:,2),1);
    else %point coordinate specification
        nodes = nodecoord(nodecoord(:,2)==BC.pos{i}{1}&nodecoord(:,3)==BC.pos{i}{2},1);
    end
    
    %BC condition cases: x,y,w,wx,wy - Displacement BCs
    %                    Nx,Ny,Nxy - Distributed Axial Load BCs (N/m)
    %                    Mx,My,Mxy - Distributed Axial Moment BCs (Nm/m)
    %                    q - Domain Transverse Load (N/m^2)
    %                    Q - Distributed Transverse Load BC (N/m)
    %all of the above could be point loads/moments
    for j = 1:size(BC.type{i},2)
        if BC.type{i}{j} == 'x'  %selected node(s)'s x displacement
            d(nodes*2-1) = BC.val{i}{j};
            d_dof(nodes*2-1) = 1;
        elseif BC.type{i}{j} == 'y'  %selected node(s)'s y displacement
            d(nodes*2) = BC.val{i}{j};
            d_dof(nodes*2) = 1;
        elseif BC.type{i}{j} == 'w'  %selected node(s)'s z displacement
            d(nodes*3-2+2*max(nodecoord(:,1))) = BC.val{i}{j};
            d_dof(nodes*3-2+2*max(nodecoord(:,1))) = 1;
        elseif strcmp(BC.type{i}{j},'wx')  %selected node(s)'s rotation about x 
            d(nodes*3-1+2*max(nodecoord(:,1))) = BC.val{i}{j};
            d_dof(nodes*3-1+2*max(nodecoord(:,1))) = 1;
        elseif strcmp(BC.type{i}{j},'wy')  %selected node(s)'s rotation about y 
            d(nodes*3+2*max(nodecoord(:,1))) = BC.val{i}{j};
            d_dof(nodes*3+2*max(nodecoord(:,1))) = 1;
        elseif strcmp(BC.type{i}{j},'Nx') %selected node(s)'s load along x-axis
            %(this should only be specified for boundary parallel to y-axis
            if length(BC.pos{i}{1}) == 1  
                BCval = [BC.val{i}{j};0]; 
                curpos = BC.pos{i};
                ax = 2; 
                shapefcn = 1;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'Ny') %selected node(s)'s load along y-axis 
            %(this should only be specified for boundary parallel to x-axis
            if length(BC.pos{i}{2}) == 1  
                BCval = [0;BC.val{i}{j}]; 
                curpos = BC.pos{i};
                ax = 1; 
                shapefcn = 1;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'Nxy') %this is for shear load parallel to y-axis (on constant x surface)
            %note the boundary surface is also parallel to y-axis
            if length(BC.pos{i}{1}) == 1  
                BCval = [0;BC.val{i}{j}]; 
                curpos = BC.pos{i};
                ax = 2; 
                shapefcn = 1;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'Nyx') %this is for shear load parallel to x-axis (on constant y surface)
            %note the boundary surface is also parallel to x-axis
            if length(BC.pos{i}{2}) == 1  
                BCval = [BC.val{i}{j};0]; 
                curpos = BC.pos{i};
                ax = 1; 
                shapefcn = 1;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'Mx') %this is for bending moment applied along the x-axis (rotate about y-axis)
            %note the boundary surface is also parallel to y-axis (constant x)
            if length(BC.pos{i}{1}) == 1  
                BCval = [BC.val{i}{j};0]; 
                curpos = BC.pos{i};
                ax = 2; 
                shapefcn = 2;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'My') %this is for bending moment applied along the y-axis (rotate about x-axis)
            %note the boundary surface is also parallel to x-axis (constant y)
            if length(BC.pos{i}{2}) == 1  
                BCval = [0;BC.val{i}{j}]; 
                curpos = BC.pos{i};
                ax = 1; 
                shapefcn = 2;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'Mxy') %this is for twisting moment applied along the y-axis (rotate about x-axis)
            %note the boundary surface is parallel to y-axis (constant x)
            if length(BC.pos{i}{1}) == 1  
                BCval = [0;BC.val{i}{j}]; 
                curpos = BC.pos{i};
                ax = 2; 
                shapefcn = 2;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j}, 'Myx') %this is for twisting moment applied along the x-axis (rotate about y-axis)
            %note the boundary surface is also parallel to x-axis (constant y)
            if length(BC.pos{i}{2}) == 1  
                BCval = [BC.val{i}{j};0]; 
                curpos = BC.pos{i};
                ax = 1; 
                shapefcn = 2;
                f = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p);
            end
        elseif strcmp(BC.type{i}{j},'Q') %this is for transverse(through-thickness shear)
            if length(BC.pos{i}{1}) == 1              
                ax = 2;
            elseif length(BC.pos{i}{2}) == 1
                ax = 1;               
            end 
            curpos = BC.pos{i};
            shapefcn = 3; 
            f = constructF_noncon(connec, nodecoord, f, curpos, BC.val{i}{j}, ax, shapefcn, p);
        elseif strcmp(BC.type{i}{j},'q')  %Domain Transverse Load (N/m^2)
            curpos = BC.pos{i};
            shapefcn = 4;
            ax = 0;
            f = constructF_noncon(connec, nodecoord, f, curpos, BC.val{i}{j}, ax, shapefcn, p);
            
        end
    end
end

fdiff = f(d_dof==0)-ktot(d_dof==0,d_dof==1)*d(d_dof==1);
d(d_dof == 0) = ktot(d_dof == 0,d_dof==0)\fdiff;
f(d_dof == 1) = ktot(d_dof==1,:)*d;
end
