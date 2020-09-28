function [J,XYZ] = fun_my_postRT2(Jout,Kt,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d)
dx=d(1); 
dy=d(2);
dz=d(3);
mytic_post_vec = tic;
J=zeros(Kt+1,3); %add fake voxel
J(1:end-1,1)=reshape(0.5*Jout(:,:,:,1)/(dy*dz),Kt,1); 
J(1:end-1,2)=reshape(0.5*Jout(:,:,:,2)/(dx*dz),Kt,1); 
J(1:end-1,3)=reshape(0.5*Jout(:,:,:,3)/(dy*dx),Kt,1); 
ind2x=Ae1x(2,:);
ind2x(ind2x==0)=Kt+1;
ind2y=Ae1y(2,:);
ind2y(ind2y==0)=Kt+1;
ind2z=Ae1z(2,:);
ind2z(ind2z==0)=Kt+1;
J(ind2x,1)=J(ind2x,1)+reshape(0.5*Jout(:,:,:,1)/(dy*dz),Kt,1);
J(ind2y,2)=J(ind2y,2)+reshape(0.5*Jout(:,:,:,2)/(dx*dz),Kt,1);
J(ind2z,3)=J(ind2z,3)+reshape(0.5*Jout(:,:,:,3)/(dy*dx),Kt,1);
mytoc_post_vec = toc(mytic_post_vec);
J(end,:)=[];
XYZ = reshape(xyz,L*M*N,3); %3D coordinates
end

