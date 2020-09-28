function [face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross] = ...
    prepareINTEGRATION(d)
%%Voxel Properties
dx = d(1); dy = d(2); dz = d(3);
%%Reference Voxel Vertices
V(1,:) = [dx/2,-dy/2,dz/2];
V(2,:) = [dx/2,-dy/2,-dz/2];
V(3,:) = [dx/2, dy/2,-dz/2];
V(4,:) = [dx/2, dy/2,dz/2];
V(5,:) = [-dx/2,-dy/2,dz/2];
V(6,:) = [-dx/2,-dy/2,-dz/2];
V(7,:) = [-dx/2,dy/2,-dz/2];
V(8,:) = [-dx/2,dy/2,dz/2];
%%Assing faces
face(:,:,1) = [V(1,:);V(2,:);V(3,:);V(4,:)]; % +x 
face(:,:,2) = [V(4,:);V(3,:);V(7,:);V(8,:)]; % +y 
face(:,:,3) = [V(1,:);V(4,:);V(8,:);V(5,:)]; % +z 
face(:,:,4) = [V(5,:);V(6,:);V(2,:);V(1,:)]; % -y
face(:,:,5) = [V(3,:);V(2,:);V(6,:);V(7,:)]; % -z
face(:,:,6) = [V(5,:);V(8,:);V(7,:);V(6,:)]; % -x
%%Face differences
face_diff = zeros(4,3,6);
for ii = 1:6
    face_diff(1,:,ii) = face(2,:,ii)-face(1,:,ii);
    face_diff(2,:,ii) = face(3,:,ii)-face(2,:,ii);
    face_diff(3,:,ii) = face(4,:,ii)-face(3,:,ii);
    face_diff(4,:,ii) = face(1,:,ii)-face(4,:,ii);
end
%%Compute r_f arbitrary point on face
%Choose a vertex
r_f = zeros(6,3);
for ii = 1:6 
    r_f(ii,:) = face(1,:,ii); %
end
%%Computation of r_e and u_e vectors for each face
r_e = zeros(4,3,6);
u_e = zeros(4,3,6);
for ii = 1:6
    for jj = 1:4     
        r_e(jj,:,ii) = face(jj,:,ii); 
    end
    u_e(1,:,ii) = (face(2,:,ii)-face(1,:,ii))/norm(face(2,:,ii)-face(1,:,ii));
    u_e(2,:,ii) = (face(3,:,ii)-face(2,:,ii))/norm(face(3,:,ii)-face(2,:,ii));
    u_e(3,:,ii) = (face(4,:,ii)-face(3,:,ii))/norm(face(4,:,ii)-face(3,:,ii));
    u_e(4,:,ii) = (face(1,:,ii)-face(4,:,ii))/norm(face(1,:,ii)-face(4,:,ii));  
end
%%Computation of triangular patch for each face 
triang_face = zeros(3,3,2,6);
for ii = 1:6
    triang_face(:,:,1,ii) = [face(1,:,ii);face(2,:,ii);face(3,:,ii)];
    triang_face(:,:,2,ii) = [face(3,:,ii);face(4,:,ii);face(1,:,ii)];
end   
%%Computation of face indices for dot product
%Which sign (+/-) and which element take in a dot product with a face
%normal. First position element index, second sign
idxf_dot = zeros(6,2);
idxf_dot(1,:) = [1,+1];
idxf_dot(2,:) = [2,+1];
idxf_dot(3,:) = [3,+1];
idxf_dot(4,:) = [2,-1];
idxf_dot(5,:) = [3,-1];
idxf_dot(6,:) = [1,-1];
%%Computation of face indices for cross product
%Which sign (+/-) and which element take in a cross product with a face
%normal. First positions element index, second poitions signs
idxf_cross = zeros(2,3,6); %as a 3D array
idxf_cross(:,:,1) = [1,3,2;0,-1,+1];
idxf_cross(:,:,2) = [3,2,1;1,0,-1];
idxf_cross(:,:,3) = [2,1,3;-1,+1,0];
idxf_cross(:,:,4) = [3,2,1;-1,0,+1];
idxf_cross(:,:,5) = [2,1,3;+1,-1,0];
idxf_cross(:,:,6) = [1,3,2;0,+1,-1];
end

