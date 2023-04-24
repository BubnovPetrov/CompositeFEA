function [KA,KB,KC,KD] = elelaminatestiff_noncon(p,x_e,y_e,ply)
%this function computes the complete [A,B;B,D] laminate stiffness
%matrix with multiplication with the shape functions and the appropriate
%derivatives

%ply: [G12,E1,E2,v12,v21,theta,h]
plynum = size(ply,1); 
%independent material constants:
G12 = ply(:,1); E1 = ply(:,2); E2 = ply(:,3); v12 = ply(:,4); v21 = ply(:,5);
theta = ply(:,6); h = ply(:,7); 
%note: v12, v21, E1, E2 are coupled: v21 = E2.*v12./E1

%material compliance matrix for plies
for i = 1:plynum
    S_o(:,:,i) = [1/E1(i),-v21(i)/E2(i),0;-v12(i)/E1(i),1/E2(i),0;0,0,1/G12(i)].*10^10;
    %D_o(:,:,i) =[E1(i)/(1-v12(i)*v21(i)),v12(i)*E2(i)/(1-v12(i)*v21(i)),0;v12(i)*E2(i)/(1-v12(i)*v21(i)),E2(i)/(1-v12(i)*v21(i)),0;0,0,G12(i)];
end

%set up variables for finding the modulus-weighted midplane 
E0 = max(E1); %reference modulus for midplane calculation
hstar = sum((E1./E0).*h);
z_hstar = 0;

%transformation of compliance matrix to principal directions for the
%laminate:
R = [1,0,0;0,1,0;0,0,2]; %Reuter's matrix 
for i = 1:plynum
    c2 = cosd(theta(i)*(180/pi))^2; 
    s2 = sind(theta(i)*(180/pi))^2;
    sc = sind(theta(i)*(180/pi))*cosd(theta(i)*(180/pi)); 
    At(:,:,i) = [c2,s2,2*sc;s2,c2,-2*sc;-sc,sc,c2-s2]; %stress rotation matrix 
    %strain rotation matrix would be RAR^-1
    S_bar = R*inv(At(:,:,i))*inv(R)*S_o(:,:,i)*At(:,:,i); %transformed compliance matrix 
    D_bar(:,:,i) = inv(S_bar./10^10); %transformed stiffness matrix 
    z_hstar = z_hstar + E1(i)/E0*h(i)*(sum(h(1:i))*2-h(i))/2;  
end

%zstar is the location of the modulus weighted midplane with respect to the
%reference coordinate system where the bottom of the laminate is defined as
%the z = 0 plane.
zstar = z_hstar/hstar; 

%forming the laminate stiffness matrix components
A = zeros(size(D_bar(:,:,1)));
B = zeros(size(D_bar(:,:,1)));
D = zeros(size(D_bar(:,:,1)));
for i = 1:plynum
    tol = max(eps(A),eps(D_bar(:,:,i)*((sum(h(1:i))-zstar)-(sum(h(1:i))-h(i)-zstar))));
    A = A+D_bar(:,:,i)*((sum(h(1:i))-zstar)-(sum(h(1:i))-h(i)-zstar));
    A(abs(A)<=100*tol) = 0;
    tol = max(eps(B),eps(D_bar(:,:,i)*((sum(h(1:i))-zstar)^2-(sum(h(1:i))-h(i)-zstar)^2)));
    B = B+0.5*D_bar(:,:,i)*((sum(h(1:i))-zstar)^2-(sum(h(1:i))-h(i)-zstar)^2);
    B(abs(B)<=100*tol) = 0;
    tol = max(eps(D),eps(D_bar(:,:,i)*((sum(h(1:i))-zstar)^3-(sum(h(1:i))-h(i)-zstar)^3)));
    D = D+(1/3)*D_bar(:,:,i)*((sum(h(1:i))-zstar)^3-(sum(h(1:i))-h(i)-zstar)^3);
    D(abs(D) <= 100*tol) = 0;
end

area = abs((x_e(3) - x_e(1)) * (y_e(3) - y_e(1))); % area of element (m^2)

% initialize a few variables
% total stiffness matrix for the element (segmented to 4 sub-matrices)
KA = zeros(8,8); 
KB = zeros(8,12);
KC = zeros(12,8);
KD = zeros(12,12);

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

J_det = area / 4; % determinant of Jacobian of transformation to (xi, eta)
dedx = 2/abs(x_e(2)-x_e(1));
dndy = 2/abs(y_e(3)-y_e(2));

%%%%%%%%%%%%%%%%%%%%% shape function computation %%%%%%%%%%%%%%%%%%%%%%%%%

Hw =  @(a,b)reshape([(a.*2.0+1.0).*(b-1.0).*(-1.0./4.0)-(a./8.0-1.0./8.0).*(b-1.0).*2.0, ...
    (a./8.0-1.0./8.0).*(b-1.0).*-2.0-(a./8.0-1.0./8.0).*(b.*2.0+1.0).*2.0, ...
    a.*(-1.0./4.0)-b./4.0-((b.*2.0+1.0).*(b-1.0))./4.0-(a.*2.0+1.0).*(a./8.0-1.0./8.0).*2.0-a.^2./4.0-b.^2./4.0+1.0./2.0, ...
    (a.*2.0-2.0).*(b-1.0).*(-1.0./4.0)-((a+1.0).*(b-1.0))./4.0,0.0, ...
    (a-1.0).^2.*(-1.0./4.0)-((a.*2.0-2.0).*(a+1.0))./4.0,0.0, ...
    (b.*2.0-2.0).*(a-1.0).*(-1.0./4.0)-((a-1.0).*(b+1.0))./4.0, ...
    (b-1.0).^2.*(-1.0./4.0)-((b.*2.0-2.0).*(b+1.0))./4.0, ...
    ((a.*2.0-1.0).*(b-1.0))./4.0+(a./8.0+1.0./8.0).*(b-1.0).*2.0, ...
    (a./8.0+1.0./8.0).*(b-1.0).*2.0+(a./8.0+1.0./8.0).*(b.*2.0+1.0).*2.0, ...
    a.*(-1.0./4.0)+b./4.0+((b.*2.0+1.0).*(b-1.0))./4.0+(a.*2.0-1.0).*(a./8.0+1.0./8.0).*2.0+a.^2./4.0+b.^2./4.0-1.0./2.0, ...
    (a.*2.0+2.0).*(b-1.0).*(-1.0./4.0)-((a-1.0).*(b-1.0))./4.0, ...
    0.0,(a+1.0).^2.*(-1.0./4.0)-((a.*2.0+2.0).*(a-1.0))./4.0, ...
    0.0,((b.*2.0-2.0).*(a+1.0))./4.0+((a+1.0).*(b+1.0))./4.0, ...
    (b-1.0).^2./4.0+((b.*2.0-2.0).*(b+1.0))./4.0, ...
    (a.*2.0-1.0).*(b+1.0).*(-1.0./4.0)-(a./8.0+1.0./8.0).*(b+1.0).*2.0, ...
    (a./8.0+1.0./8.0).*(b+1.0).*-2.0-(a./8.0+1.0./8.0).*(b.*2.0-1.0).*2.0, ...
    a./4.0+b./4.0-((b.*2.0-1.0).*(b+1.0))./4.0-(a.*2.0-1.0).*(a./8.0+1.0./8.0).*2.0-a.^2./4.0-b.^2./4.0+1.0./2.0, ...
    ((a.*2.0+2.0).*(b+1.0))./4.0+((a-1.0).*(b+1.0))./4.0, ...
    0.0,(a+1.0).^2./4.0+((a.*2.0+2.0).*(a-1.0))./4.0,0.0, ...
    ((b.*2.0+2.0).*(a+1.0))./4.0+((a+1.0).*(b-1.0))./4.0, ...
    (b+1.0).^2./4.0+((b.*2.0+2.0).*(b-1.0))./4.0, ...
    ((a.*2.0+1.0).*(b+1.0))./4.0+(a./8.0-1.0./8.0).*(b+1.0).*2.0, ...
    (a./8.0-1.0./8.0).*(b+1.0).*2.0+(a./8.0-1.0./8.0).*(b.*2.0-1.0).*2.0, ...
    a./4.0-b./4.0+((b.*2.0-1.0).*(b+1.0))./4.0+(a.*2.0+1.0).*(a./8.0-1.0./8.0).*2.0+a.^2./4.0+b.^2./4.0-1.0./2.0, ...
    ((a.*2.0-2.0).*(b+1.0))./4.0+((a+1.0).*(b+1.0))./4.0,0.0, ...
    (a-1.0).^2./4.0+((a.*2.0-2.0).*(a+1.0))./4.0,0.0,(b.*2.0+2.0).*(a-1.0).*(-1.0./4.0)-((a-1.0).*(b-1.0))./4.0, ...
    (b+1.0).^2.*(-1.0./4.0)-((b.*2.0+2.0).*(b-1.0))./4.0],[3,12]);

% assemble stiffness matrix
for i = 1:p %loop over eta gauss points
     for j = 1:p %loop over nu gauss points
         % current values of x and y at Gauss points
         x = 0.25*[(1-eta(i))*(1-eta(j)),(1+eta(i))*(1-eta(j)),(1+eta(i))*(1+eta(j)),(1-eta(i))*(1+eta(j))]*x_e(:);
         y = 0.25*[(1-eta(i))*(1-eta(j)),(1+eta(i))*(1-eta(j)),(1+eta(i))*(1+eta(j)),(1-eta(i))*(1+eta(j))]*y_e(:);
         % value of H at current Gauss point H is a 3 by 8 matrix
         Hsub = (1 / area) * [(y - y_e(4)), 0, -(y - y_e(4)), 0, (y - y_e(1)), 0 , - (y - y_e(1)), 0;
             0, (x - x_e(2)), 0, -(x - x_e(1)), 0, (x - x_e(1)), 0, -(x - x_e(2)); 
             (x - x_e(2)), (y - y_e(4)), -(x - x_e(1)), -(y - y_e(4)), (x - x_e(1)), (y - y_e(1)), -(x - x_e(2)), -(y - y_e(1))];  
         
         % value of Hw at current Gauss point is a 3 by 12 matrix
         Hwsub = Hw(eta(i),eta(j));
         Hwsub(1,:) = Hwsub(1,:)*dedx^2; 
         Hwsub(2,:) = Hwsub(2,:)*dndy^2; 
         Hwsub(3,:) = Hwsub(3,:)*dndy*dedx*2; 
         
        

         % contribution to stiffness matrix from current Gauss point
         tol = max(eps(KA),eps(double(w(i) * w(j) * J_det * Hsub' * A * Hsub)));
         KA = KA + double(w(i) * w(j) * J_det * Hsub' * A * Hsub);   %12 by 12
         KA(abs(KA) <= 100*tol) = 0;         

         tol = max(eps(KB),eps(double(w(i) * w(j) * J_det * Hsub' * B * Hwsub)));
         KB = KB - double(w(i) * w(j) * J_det * Hsub' * B * Hwsub);   %8 by 12
         KB(abs(KB) <= 100*tol) = 0;
         
         tol = max(eps(KC),eps(double(w(i) * w(j) * J_det * Hwsub' * B * Hsub)));
         KC = KC - double(w(i) * w(j) * J_det * Hwsub' * B * Hsub);   %12 by 8
         KC(abs(KC) <= 100*tol) = 0;
         
         tol = max(eps(KD),eps(double(w(i) * w(j) * J_det * Hwsub' * D * Hwsub)));
         KD = KD + double(w(i) * w(j) * J_det * Hwsub' * D * Hwsub);   %12 by 12
         KD(abs(KD) <= 100*tol) = 0;
     end
end

end  



