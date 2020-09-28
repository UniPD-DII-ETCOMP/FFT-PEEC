function [y] = Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,P_target,w)
%%Initializationn
Np = size(P_target,1); %number of integration points on volume
GG = zeros(6,Np);
for kk = 1:Np 
    %Subs for calculation of we
    subsw = face - P_target(kk,:);
    %%Subs for calculation of D
    subsD = triang_face - P_target(kk,:);
    for ii = 1:6
        %%Compute we eq. (18)
        epsilon = norm(face_diff(1,:,ii))/(norm(subsw(2,:,ii))+norm(subsw(1,:,ii)));
        w_e(1) = log((1+epsilon)/(1-epsilon)); 
        epsilon = norm(face_diff(2,:,ii))/(norm(subsw(3,:,ii))+norm(subsw(2,:,ii)));
        w_e(2) = log((1+epsilon)/(1-epsilon));  
        epsilon = norm(face_diff(3,:,ii))/(norm(subsw(4,:,ii))+norm(subsw(3,:,ii)));
        w_e(3) = log((1+epsilon)/(1-epsilon));
        epsilon = norm(face_diff(4,:,ii))/(norm(subsw(1,:,ii))+norm(subsw(4,:,ii)));
        w_e(4) = log((1+epsilon)/(1-epsilon));
        %%Compute D eq. (21)
        omega = zeros(2,1);
        for jj = 1:2   
            D1 = norm(subsD(1,:,jj,ii))*norm(subsD(2,:,jj,ii))*norm(subsD(3,:,jj,ii));
            D2 = norm(subsD(3,:,jj,ii))*sum(subsD(1,:,jj,ii).*subsD(2,:,jj,ii));
            D3 = norm(subsD(2,:,jj,ii))*sum(subsD(1,:,jj,ii).*subsD(3,:,jj,ii));
            D4 = norm(subsD(1,:,jj,ii))*sum(subsD(2,:,jj,ii).*subsD(3,:,jj,ii));
            D = D1 + D2 + D3 + D4;
            crv = cross(subsD(2,:,jj,ii),subsD(3,:,jj,ii));
            dt = dot(subsD(1,:,jj,ii),crv);
            %dt = my_dot(triang_face(1,:,jj,ii)-P_target(kk,:),crv);
            omega(jj) = 2*atan2(dt,D);      
        end  
        sub = r_e(:,:,ii)-repmat(P_target(kk,:),4,1);
        crn = sub(:,idxf_cross(1,:,ii)).*idxf_cross(2,:,ii);
        Omega = sum(omega);%omega(1) + omega(2);
        dt = dot(crn,u_e(:,:,ii),2); %a
        sub = r_f(ii,:) - P_target(kk,:);
        dtn = sub(idxf_dot(ii,1))*idxf_dot(ii,2); %
        W_f =  dot(dt,w_e) - (dtn*Omega);
        GG(ii,kk) = dtn*W_f;
    end    
end
I = 0.5*sum(GG,1); %sum rows
%% Numerical Integration
y = sum(w.*I');
%% Add Constants
y = y/(4*pi);
end

