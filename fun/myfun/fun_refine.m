function [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz_ori4D,smeshx_ori,smeshy_ori,smeshz_ori,Nmat_ori,L_ori,M_ori,N_ori,who,mymod)
xyz_ori=reshape(xyz_ori4D,L_ori*M_ori*N_ori,3);
% figure
% plot3(xyz_ori(:,1),xyz_ori(:,2),xyz_ori(:,3),'.')
% axis equal
xyz=zeros(L_ori*M_ori*N_ori*2,3);
if who==1
    xyz(1:L_ori*M_ori*N_ori,1)=xyz_ori(:,1)-smeshx_ori/4;
    xyz(1:L_ori*M_ori*N_ori,2)=xyz_ori(:,2);
    xyz(1:L_ori*M_ori*N_ori,3)=xyz_ori(:,3);
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,1)=xyz_ori(:,1)+smeshx_ori/4;
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,2)=xyz_ori(:,2);
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,3)=xyz_ori(:,3);
    smeshx=smeshx_ori/2;
    smeshy=smeshy_ori;
    smeshz=smeshz_ori;
    L=L_ori*2;
    M=M_ori;
    N=N_ori;
elseif who==2
    xyz(1:L_ori*M_ori*N_ori,1)=xyz_ori(:,1);
    xyz(1:L_ori*M_ori*N_ori,2)=xyz_ori(:,2)-smeshy_ori/4;
    xyz(1:L_ori*M_ori*N_ori,3)=xyz_ori(:,3);
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,1)=xyz_ori(:,1);
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,2)=xyz_ori(:,2)+smeshy_ori/4;
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,3)=xyz_ori(:,3);
    smeshx=smeshx_ori;
    smeshy=smeshy_ori/2;
    smeshz=smeshz_ori;
    L=L_ori;
    M=M_ori*2;
    N=N_ori;
elseif who==3
    xyz(1:L_ori*M_ori*N_ori,1)=xyz_ori(:,1);
    xyz(1:L_ori*M_ori*N_ori,2)=xyz_ori(:,2);
    xyz(1:L_ori*M_ori*N_ori,3)=xyz_ori(:,3)-smeshz_ori/4;
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,1)=xyz_ori(:,1);
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,2)=xyz_ori(:,2);
    xyz(L_ori*M_ori*N_ori+1:L_ori*M_ori*N_ori*2,3)=xyz_ori(:,3)+smeshz_ori/4;
    smeshx=smeshx_ori;
    smeshy=smeshy_ori;
    smeshz=smeshz_ori/2;
    L=L_ori;
    M=M_ori;
    N=N_ori*2;    
end
xmin=min(xyz(:,1));
ymin=min(xyz(:,2));
zmin=min(xyz(:,3));
xmax=max(xyz(:,1));
ymax=max(xyz(:,2));
zmax=max(xyz(:,3));
% reordering
IDtriplet=round...
            ([(xyz(:,1)-xmin+smeshx)./smeshx,...
              (xyz(:,2)-ymin+smeshy)./smeshy,...
              (xyz(:,3)-zmin+smeshz)./smeshz]);
Id_lin=IDtriplet(:,1)+(IDtriplet(:,2)-1)*L+(IDtriplet(:,3)-1)*L*M; %   
xyz(Id_lin,:)=xyz;
% figure
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.')
% axis equal
% 
Nmat=Nmat_ori;
for ii = 1:Nmat_ori
   if strcmp(Ind(ii).tag,'air')
       Ind(ii).ind=[Id_lin(Ind(ii).ind);...
                    Id_lin(L_ori*M_ori*N_ori+Ind(ii).ind)];
   elseif strcmp(Ind(ii).tag,'cond')
       Ind(ii).ind=[Id_lin(Ind(ii).ind);...
                    Id_lin(L_ori*M_ori*N_ori+Ind(ii).ind)];  
   elseif strcmp(Ind(ii).tag,'diel')
       Ind(ii).ind=[Id_lin(Ind(ii).ind);...
                    Id_lin(L_ori*M_ori*N_ori+Ind(ii).ind)];                 
   elseif strcmp(Ind(ii).tag,'mag')
       Ind(ii).ind=[Id_lin(Ind(ii).ind);...
                    Id_lin(L_ori*M_ori*N_ori+Ind(ii).ind)];                   
   elseif strcmp(Ind(ii).tag,'terminal') || strcmp(Ind(ii).tag,'cur')
       if mymod==1 % ports only in one of the two splitted voxels
           tmp=(Ind(ii).ind);
           Ind(ii).ind=Id_lin(tmp);
           Ind(Nmat+1).ind = Id_lin(L_ori*M_ori*N_ori+tmp);  
           Ind(Nmat+1).tag = 'cond';
           Ind(Nmat+1).rho = Ind(ii).rho;
           Nmat=Nmat+1;
       elseif mymod==2  % ports in both the two splitted voxels
           Ind(ii).ind=[Id_lin(Ind(ii).ind);...
                    Id_lin(L_ori*M_ori*N_ori+Ind(ii).ind)];             
       end
   end
end
xyz=reshape(xyz,L,M,N,3);
nVoxel = L*M*N;
end

