function [fout] = constructF_noncon(connec, nodecoord, f, curpos, BCval, ax, shapefcn, p)
%for non-conforming case
%input var: connec, nodecoord, nodes, f, curpos = BC.pos{i}, ax,
%shapefcn,p, BCval

%connec, nodecoord, nodes are taken from solve as it is 

%f is the last updated force vector

%BCval is a syms or value vector/scalar 

%curpos = BC.pos{i} from solve, it's a cell vector 

%p: Gauss quadrature order

%shapefcn: associated shape function 
%shapefcn = 1: in-plane stress loading(Nx,Ny,Nxy)
%shapefcn = 2: moments(Mx,My,Mxy)
%shapefcn = 3 and 4: transverse loading (Q,q)

%ax: integration boundary 
%ax = 0: area of integration 
%ax = 1: boundary parallel to x-axis (constant y)
%ax = 2: boundary parallel to y-axis (constant x)

%find all the relevant elements 
elem = []; 
ispoint = 0;
%for i = 1:length(nodes) %loop over all the relevant BC nodes

%assume the applied BC is either a point load/moment, or the length scale
%of the BC is at least the size of one element's dimension 
for i = 1:size(connec,1) %loop over all elements to find relevant ones
    %val = nodes(i); %current node
    %find all elements attached to the node
    curelem = connec(i,:); 
    flag = 0;
    if sum(curelem(1)==elem)==0 %make sure element is not a duplicate
        if length(curpos{1}) == 1 && length(curpos{2})==1 %point load case 
            ispoint = 1;
            if curpos{1} >= nodecoord(curelem(2),2) && ... 
                    curpos{1} <= nodecoord(curelem(3),2) && ... 
                    curpos{2} >= nodecoord(curelem(2),3) && ...
                    curpos{2} <= nodecoord(curelem(4),3) 
                %the point load application is with in the current element
            else
                %the point load is not in the current element 
                flag = 1; 
            end
        else
   
            if ax == 0 %distributed area loading 
                %if at least 1 node connected to the element are within BC
                %application area, then leave flag to be 0 
                if nodecoord(curelem(2),2) <= curpos{1}(end) && nodecoord(curelem(2),3) <= curpos{2}(end) ...
                    && nodecoord(curelem(2),2) >= curpos{1}(1) && nodecoord(curelem(2),3) >= curpos{2}(1)
                elseif nodecoord(curelem(3),2) <= curpos{1}(end) && nodecoord(curelem(3),3) <= curpos{2}(end) ...
                    && nodecoord(curelem(3),2) >= curpos{1}(1) && nodecoord(curelem(3),3) >= curpos{2}(1)
                elseif nodecoord(curelem(3),2) <= curpos{1}(end) && nodecoord(curelem(3),3) <= curpos{2}(end) ...
                    && nodecoord(curelem(3),2) >= curpos{1}(1) && nodecoord(curelem(3),3) >= curpos{2}(1)
                elseif nodecoord(curelem(4),2) <= curpos{1}(end) && nodecoord(curelem(4),3) <= curpos{2}(end) ...
                    && nodecoord(curelem(4),2) >= curpos{1}(1) && nodecoord(curelem(4),3) >= curpos{2}(1)
                else
                    %if none of nodes connected to the element is inside of BC
                    %boundary, then change flag to be 1
                    flag = 1;
                end
            elseif ax == 1 %distributed boundary line loading along x (constant y surface)
                if nodecoord(curelem(2),3) == curpos{2} || nodecoord(curelem(4),3) == curpos{2}
                    if nodecoord(curelem(2),2) <= curpos{1}(end) && nodecoord(curelem(2),2) >= curpos{1}(1)
                        %if the any of the nodes connected to the element is within BC
                        %x-axis boundary, then leave flag to be 0
                    elseif nodecoord(curelem(3),2) <= curpos{1}(end) && nodecoord(curelem(3),2) >= curpos{1}(1)
                    else
                        flag = 1;
                    end
                else
                    %otherwise, change flag to be 1
                    flag = 1;
                end
            elseif ax == 2 %distributed boundary line loading along y (constant x surface)
                if nodecoord(curelem(2),2) == curpos{1} || nodecoord(curelem(4),2) == curpos{1}
                    if nodecoord(curelem(2),3) <= curpos{2}(end) && nodecoord(curelem(2),3) >= curpos{2}(1)
                    %if the all nodes connected to the element is within BC
                    %boundary, then leave flag to be 0
                    elseif nodecoord(curelem(4),3) <= curpos{2}(end) && nodecoord(curelem(4),3) >= curpos{2}(1)
                    else
                        flag = 1;
                    end
                else
                    %if the any nodes connected to the element is outside of BC
                    %boundary, then change flag to be 1
                    flag = 1;
                end
            end
            %end
        end
        %if all nodes of the element is within BC boundary, then this
        %element will be apart of the Gauss quadrature computation
        if flag == 0
            elem = [elem;curelem(1)];
        end    
    end
    %end
end

% set up Gauss quadrature
if p == 1
    eta = 0; %quadrature points 
    w = 2; %quadrature weights
elseif p == 2
    eta = [-1/sqrt(3),1/sqrt(3)];
    w = [1,1];
elseif p == 3
    eta = [0,-sqrt(3/5),sqrt(3/5)];
    w = [8/9,5/9,5/9];
elseif p == 4
    eta = [-sqrt(3/7 - (2/7)*sqrt(6/5)),sqrt(3/7 - (2/7)*sqrt(6/5)),-sqrt(3/7 + (2/7)*sqrt(6/5)),sqrt(3/7 + (2/7)*sqrt(6/5))];
	w = [(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36];
end

%%%%%%%%%%%%%%%%%%%%% shape function computation %%%%%%%%%%%%%%%%%%%%%%%%%
%shapefcn = 1: in-plane stress loading(Nx,Ny,Nxy)
%shapefcn = 2: moments(Mx,My,Mxy)
%shapefcn = 3 and 4: transverse loading (Q,q)

if shapefcn == 2
    B = @(a,b)reshape([(b-1.0).*(a+b+a.^2+b.^2-2.0).*(-1.0./8.0)-(a.*2.0+1.0).*(a./8.0-1.0./8.0).*(b-1.0), ...
        -(a./8.0-1.0./8.0).*(a+b+a.^2+b.^2-2.0)-(a./8.0-1.0./8.0).*(b.*2.0+1.0).*(b-1.0),(a-1.0).^2.*(b-1.0).*(-1.0./8.0)-((a.*2.0-2.0).*(a+1.0).*(b-1.0))./8.0, ...
        (a-1.0).^2.*(a+1.0).*(-1.0./8.0),(b-1.0).^2.*(b+1.0).*(-1.0./8.0),(a-1.0).*(b-1.0).^2.*(-1.0./8.0)-((b.*2.0-2.0).*(a-1.0).*(b+1.0))./8.0, ...
        ((b-1.0).*(-a+b+a.^2+b.^2-2.0))./8.0+(a.*2.0-1.0).*(a./8.0+1.0./8.0).*(b-1.0), ...
        (a./8.0+1.0./8.0).*(-a+b+a.^2+b.^2-2.0)+(a./8.0+1.0./8.0).*(b.*2.0+1.0).*(b-1.0), ...
        (a+1.0).^2.*(b-1.0).*(-1.0./8.0)-((a.*2.0+2.0).*(a-1.0).*(b-1.0))./8.0, ...
        (a-1.0).*(a+1.0).^2.*(-1.0./8.0),((b-1.0).^2.*(b+1.0))./8.0, ...
        ((a+1.0).*(b-1.0).^2)./8.0+((b.*2.0-2.0).*(a+1.0).*(b+1.0))./8.0, ...
        ((b+1.0).*(a+b-a.^2-b.^2+2.0))./8.0-(a.*2.0-1.0).*(a./8.0+1.0./8.0).*(b+1.0), ...
        (a./8.0+1.0./8.0).*(a+b-a.^2-b.^2+2.0)-(a./8.0+1.0./8.0).*(b.*2.0-1.0).*(b+1.0), ...
        ((a+1.0).^2.*(b+1.0))./8.0+((a.*2.0+2.0).*(a-1.0).*(b+1.0))./8.0, ...
        ((a-1.0).*(a+1.0).^2)./8.0,((b-1.0).*(b+1.0).^2)./8.0, ...
        ((a+1.0).*(b+1.0).^2)./8.0+((b.*2.0+2.0).*(a+1.0).*(b-1.0))./8.0, ...
        ((b+1.0).*(a-b+a.^2+b.^2-2.0))./8.0+(a.*2.0+1.0).*(a./8.0-1.0./8.0).*(b+1.0), ...
        (a./8.0-1.0./8.0).*(a-b+a.^2+b.^2-2.0)+(a./8.0-1.0./8.0).*(b.*2.0-1.0).*(b+1.0), ...
        ((a-1.0).^2.*(b+1.0))./8.0+((a.*2.0-2.0).*(a+1.0).*(b+1.0))./8.0,((a-1.0).^2.*(a+1.0))./8.0, ...
        (b-1.0).*(b+1.0).^2.*(-1.0./8.0),(a-1.0).*(b+1.0).^2.*(-1.0./8.0)-((b.*2.0+2.0).*(a-1.0).*(b-1.0))./8.0],[2,12]);

elseif shapefcn == 3 || shapefcn == 4
    B = @(a,b)[-(a./8.0-1.0./8.0).*(b-1.0).*(a+b+a.^2+b.^2-2.0), ...
        (a-1.0).^2.*(a+1.0).*(b-1.0).*(-1.0./8.0),(a-1.0).*(b-1.0).^2.*(b+1.0).*(-1.0./8.0), ...
        (a./8.0+1.0./8.0).*(b-1.0).*(-a+b+a.^2+b.^2-2.0),(a-1.0).*(a+1.0).^2.*(b-1.0).*(-1.0./8.0), ...
        ((a+1.0).*(b-1.0).^2.*(b+1.0))./8.0,(a./8.0+1.0./8.0).*(b+1.0).*(a+b-a.^2-b.^2+2.0), ...
        ((a-1.0).*(a+1.0).^2.*(b+1.0))./8.0,((a+1.0).*(b-1.0).*(b+1.0).^2)./8.0, ...
        (a./8.0-1.0./8.0).*(b+1.0).*(a-b+a.^2+b.^2-2.0),((a-1.0).^2.*(a+1.0).*(b+1.0))./8.0, ...
        (a-1.0).*(b-1.0).*(b+1.0).^2.*(-1.0./8.0)];

elseif shapefcn == 1
    B = @(a,b)reshape([(a./4.0-1.0./4.0).*(b-1.0),0.0,0.0,(a./4.0-1.0./4.0).*(b-1.0), ...
        -(a./4.0+1.0./4.0).*(b-1.0),0.0,0.0,-(a./4.0+1.0./4.0).*(b-1.0),(a./4.0+1.0./4.0).*(b+1.0),0.0,0.0, ...
        (a./4.0+1.0./4.0).*(b+1.0),-(a./4.0-1.0./4.0).*(b+1.0),0.0,0.0,-(a./4.0-1.0./4.0).*(b+1.0)],[2,8]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = f; 
if ispoint == 1 %case of point load/moment
    split = BCval./(length(elem)); %redistribute the original value to all the attached elements 
end
%Populate fout
syms x y
for i = 1:length(elem) %loop over all relevant elements
    x_e = nodecoord(connec(elem(i),2:5),2); 
    y_e = nodecoord(connec(elem(i),2:5),3); 
    dedx = 2/abs(x_e(2)-x_e(1));
    dndy = 2/abs(y_e(3)-y_e(2));
    rowpos1 = [connec(elem(i),2)*2-1,connec(elem(i),2)*2,connec(elem(i),3)*2-1,connec(elem(i),3)*2,connec(elem(i),4)*2-1,connec(elem(i),4)*2,connec(elem(i),5)*2-1,connec(elem(i),5)*2];
    rowpos2 = [connec(elem(i),2)*3-2,connec(elem(i),2)*3-1,connec(elem(i),2)*3,connec(elem(i),3)*3-2,connec(elem(i),3)*3-1,connec(elem(i),3)*3, ... 
            connec(elem(i),4)*3-2,connec(elem(i),4)*3-1,connec(elem(i),4)*3,connec(elem(i),5)*3-2,connec(elem(i),5)*3-1,connec(elem(i),5)*3]+ 2*max(nodecoord(:,1));
    if ispoint == 0 %not a point load 
        if ax == 0 %the case of distributed q load 
            ftemp = zeros(12,1);
            if isa(BCval,'sym')
                for k = 1:p %loop over eta gauss points
                     for j = 1:p %loop over nu gauss points
                         % current values of x and y at Gauss points
                         xc = 0.25*[(1-eta(k))*(1-eta(j)),(1+eta(k))*(1-eta(j)),(1+eta(k))*(1+eta(j)),(1-eta(k))*(1+eta(j))]*x_e(:);
                         yc = 0.25*[(1-eta(k))*(1-eta(j)),(1+eta(k))*(1-eta(j)),(1+eta(k))*(1+eta(j)),(1-eta(k))*(1+eta(j))]*y_e(:);
                         J_det = (1/dedx)*(1/dndy);
                         % contribution to fout vector from current Gauss point
                         if xc < curpos{1}(1) || xc > curpos{1}(end) || yc < curpos{2}(1) || yc > curpos{2}(end)
                         else
                            Bsub = B(eta(k),eta(j));
                            ftemp = ftemp - double(w(k) * w(j) * J_det * Bsub' * subs(BCval,{x,y},{xc,yc}));
                         end
                     end
                end
            else
                for k = 1:p %loop over eta gauss points
                     for j = 1:p %loop over nu gauss points
                         xc = 0.25*[(1-eta(k))*(1-eta(j)),(1+eta(k))*(1-eta(j)),(1+eta(k))*(1+eta(j)),(1-eta(k))*(1+eta(j))]*x_e(:);
                         yc = 0.25*[(1-eta(k))*(1-eta(j)),(1+eta(k))*(1-eta(j)),(1+eta(k))*(1+eta(j)),(1-eta(k))*(1+eta(j))]*y_e(:);
                         J_det = (1/dedx)*(1/dndy);
                         %J_det = 1;
                         % contribution to fout vector from current Gauss point
                         if xc < curpos{1}(1) || xc > curpos{1}(end) || yc < curpos{2}(1) || yc > curpos{2}(end)
                         else
                            Bsub = B(eta(k),eta(j));
                            ftemp = ftemp + (1*BCval)*w(k) * w(j) * J_det * Bsub';   %12 by 1
                         end
                     end
                end
            end
            fout(rowpos2) = fout(rowpos2)+ftemp;
        elseif ax==1 %boundary parallel to x-axis (constant y surface)
            if shapefcn == 1
                ftemp = zeros(8,1);
            elseif shapefcn == 3 || shapefcn == 2
                ftemp = zeros(12,1);
            end
            yc = curpos{2}; %y-value is constant
            nu = (2*curpos{2}-y_e(1)-y_e(4))/(y_e(4)-y_e(1));
            % current values of x and y at Gauss points
            J_det = (1/dedx);
            if isa(BCval,'sym')
                for k = 1:p %loop over eta gauss points
                     % contribution to fout vector from current Gauss point
                     xc = 0.25*[1-eta(k),1+eta(k),1+eta(k),1-eta(k)]*x_e(:);
                     if xc < curpos{1}(1) || xc > curpos{1}(end)
                     else
                         if shapefcn == 2
                             %Bsub = subs(B,{a,b},{eta(k),nu});
                             Bsub = B(eta(k),nu);
                             Bsub(1,:) = Bsub(1,:)*dedx;
                             Bsub(2,:) = Bsub(2,:)*dndy;
                         else
                             %Bsub = subs(B,{a,b},{eta(k),nu});
                             Bsub = B(eta(k),nu);
                         end
                         ftemp = ftemp + double(w(k) * J_det * Bsub' * subs(BCval,{x,y},{xc,yc})); 
                     end
                end
            else
                for k = 1:p %loop over eta gauss points
                     % contribution to fout vector from current Gauss point
                     xc = 0.25*[1-eta(k),1+eta(k),1+eta(k),1-eta(k)]*x_e(:);
                     if xc < curpos{1}(1) || xc > curpos{1}(end)
                     else
                         if shapefcn == 2
                             Bsub = B(eta(k),nu);
                             Bsub(1,:) = Bsub(1,:)*dedx;
                             Bsub(2,:) = Bsub(2,:)*dndy;
                         else
                             Bsub = B(eta(k),nu);
                         end
                         ftemp = ftemp + double(w(k)* J_det * Bsub')*BCval;
                     end
                end
            end
            if shapefcn == 1
                fout(rowpos1) = fout(rowpos1)+ftemp;
            elseif shapefcn == 3 || shapefcn == 2
                fout(rowpos2) = fout(rowpos2)+ftemp;
            end
        elseif ax == 2 %boundary parallel to y-axis 
            if shapefcn == 1
                ftemp = zeros(8,1);
            elseif shapefcn == 3 || shapefcn == 2
                ftemp = zeros(12,1);
            end
            xc = curpos{1}; %x-value is constant
            nu = (2*curpos{1}-x_e(1)-x_e(2))/(x_e(2)-x_e(1));
            J_det = (1/dndy);
            if isa(BCval,'sym')
                for k = 1:p %loop over eta gauss points
                     % current values of x and y at Gauss points
                     yc = 0.25*[1-eta(k),1+eta(k),1+eta(k),1-eta(k)]*y_e(:);
                     if yc < curpos{2}(1) || yc > curpos{2}(end)
                     else
                         if shapefcn == 2
                             Bsub = B(nu,eta(k));
                             Bsub(1,:) = Bsub(1,:)*dedx;
                             Bsub(2,:) = Bsub(2,:)*dndy;
                         else
                             Bsub = B(nu,eta(k));
                         end
                         % contribution to fout vector from current Gauss point
                         ftemp = ftemp + double(w(k) * J_det * Bsub' * subs(BCval,{x,y},{xc,yc}));   
                     end
                end
            else
                for k = 1:p %loop over eta gauss points
                    yc = 0.25*[1-eta(k),1+eta(k),1+eta(k),1-eta(k)]*y_e(:);
                    if yc < curpos{2}(1) || yc > curpos{2}(end)
                    else
                         if shapefcn == 2
                             Bsub = B(nu,eta(k));
                             Bsub(1,:) = Bsub(1,:)*dedx;
                             Bsub(2,:) = Bsub(2,:)*dndy;
                         else
                             Bsub = B(nu,eta(k));
                         end

                         % contribution to fout vector from current Gauss point
                         ftemp = ftemp + double(w(k)* J_det * Bsub')*BCval; 
                    end
                end
            end
            if shapefcn == 1
                fout(rowpos1) = fout(rowpos1)+ftemp;
            elseif shapefcn == 3 || shapefcn == 2
                fout(rowpos2) = fout(rowpos2)+ftemp;
            end
        end
    elseif ispoint == 1 %point load/moment
        if shapefcn == 1
            ftemp = zeros(8,1);
        elseif shapefcn == 3 || shapefcn == 2 || shapefcn == 4
            ftemp = zeros(12,1);
        end
        yc = curpos{2}; %y-value is constant
        nuc = (2*yc-y_e(1)-y_e(4))/(y_e(4)-y_e(1));  
        xc = curpos{1}; %x-value is constant
        etac = (2*xc-x_e(1)-x_e(2))/(x_e(2)-x_e(1));
        Bsub = B(etac,nuc);
        if shapefcn == 2 %case of point moment M
            Bsub(1,:) = Bsub(1,:)*dedx;
            Bsub(2,:) = Bsub(2,:)*dndy;
        end
        if shapefcn == 4  %case of point transverse load q 
            ftemp = ftemp - double(Bsub' * split); 
        elseif shapefcn == 3 || shapefcn == 1 || shapefcn == 2
            %case of point transverse shear Q, point in-plane load N, point moment M
            ftemp = ftemp + double(Bsub' * split); 
        end
        if shapefcn == 1
            fout(rowpos1) = fout(rowpos1)+ftemp;
        elseif shapefcn == 3 || shapefcn == 2 || shapefcn == 4
            fout(rowpos2) = fout(rowpos2)+ftemp;
        end
    end
end

end