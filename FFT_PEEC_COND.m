%% FFT-PEEC: CONDUCTORS
close all
clear global
clear
clc
restoredefaultpath
warning on
format short
%%
%% BEGIN USER SETTINGS
%%
%% Directory
name_dir='test6';
%% Frequency
freq = 5e6; %[Hz]
%% Selections
plot_vectorsJ_flag = 1; %quiver plot of real and imag of J
plot_potential_flag = 1; %color plot of phi real and imag
paraview_export_flag = 1; % export to paraviw
refine.flag = 0; refine.x=1; refine.y=1; refine.z=1; % refine
capacitive_effects_flag = 1;  % capacitive effects flag 
retardation_flag = 0;         % retarded potential flag
Integration_flag = 'NumAn'; %'NumAn'; 'NumNum' (Integration: NumericalNumerical or AnalyticalNumerical)
ext_field_flag = 0; % exernal field
% below you can write the external electric field as a function of x,y,z
% and omega. Active only if ext_field_flag=1
Ex_ext = @(x,y,z,omega) -1j*omega*y; Ey_ext = @(x,y,z,omega) 1j*omega*x; Ez_ext = @(x,y,z,omega) 0*z; % external field 
%% Solver parameters
tol = 1e-6;
inner_it = 90;
outer_it = 1;
%%
%% END USER SETTINGS
%%
%% Add Path
dad = pwd;
cd('fun'); addpath(genpath(pwd)); cd(dad)
cd('fortran'); addpath(pwd); cd(dad)
cd('data'); cd(name_dir); load('data.mat'); 
fileList = dir('*.stl');
figure
hold on
xmin=[];xmax=[];ymin=[];ymax=[];zmin=[];zmax=[];ccolor=distinguishable_colors(size(fileList,1));
for ii = 1:size(fileList,1)
    [stlcoords] = READ_stl(fileList(ii).name);
    xco = squeeze( stlcoords(:,1,:) )';
    yco = squeeze( stlcoords(:,2,:) )';
    zco = squeeze( stlcoords(:,3,:) )';
    [hpat] = patch(xco,yco,zco,ccolor(ii,:));
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(3)
    title('stl (original, not scaled)')
    drawnow
end
try
   disp('-------------------------------------------------------------------')
   disp(' ') 
   type('comments.txt') 
   disp(' ')
   disp('PRESS ANY BUTTON TO CONTINUE')
   disp(' ')
   pause
   disp(' ')
   disp('-------------------------------------------------------------------')
end
cd(dad)
modelname = name_dir;
%% refine
if refine.flag
    warning('refine on')
    mymod=1;
    for ii = 1:refine.x
        [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz,smeshx,smeshy,smeshz,Nmat,L,M,N,1,mymod);
    end
    for ii = 1:refine.y
        [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz,smeshx,smeshy,smeshz,Nmat,L,M,N,2,mymod);
    end   
    for ii = 1:refine.z
        [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz,smeshx,smeshy,smeshz,Nmat,L,M,N,3,mymod);
    end      
end
%% EM constants
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*freq;
%% extract data information 
rhoVoxel=zeros(nVoxel,1);
idxV=[]; rhomin=Inf; ind_c=[]; val_c=[]; k=1;
for ii = 1:Nmat
    Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
    if strcmp(Ind(ii).tag,'air') || strcmp(Ind(ii).tag,'mag') || strcmp(Ind(ii).tag,'diel')
        % nothing to do here (?)
    elseif strcmp(Ind(ii).tag,'cond')
        idxV=[idxV;Ind(ii).ind];  
        rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;  
        rhomin=min([rhomin,Ind(ii).rho]);
    elseif strcmp(Ind(ii).tag,'terminal')
        idxV=[idxV;Ind(ii).ind]; 
        rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;
        rhomin=min([rhomin,Ind(ii).rho]);
        ind_c(k)=Ind(ii).ind(1);
        val_c(k)=Ind(ii).cur;
        k=k+1;
    end
end
idxV=unique(idxV);
del = sqrt(2*rhomin/omega/mu); %skin effect: penetration depth
%% Grid Definition
disp('----DOMAIN--------------------------------')
%%Grid resolution
disp([' Number of voxels in x direction: ', num2str(L)])
disp([' Number of voxels in y direction: ', num2str(M)])
disp([' Number of voxels in z direction: ', num2str(N)])
disp(' Resolution:')
dx = smeshx; dy = smeshy; dz = smeshz;
disp([' dx = ',num2str(dx),' m']); disp([' dy = ',num2str(dy),' m']); disp([' dz = ',num2str(dz),' m'])
d = [dx dy dz]; 
Kt = nVoxel; %total number of voxels
K = length(idxV); %number of non-empty voxels
%% Set Material Properties
rho_eV=reshape(rhoVoxel,L,M,N); % 
clear rhoVoxel
%%
disp([' Total number of voxels: ', num2str(Kt)])
disp([' Number of non-empty voxels: ', num2str(K)])
disp(' ')
%% Incidence Matix A
disp('----COMPUTING INCIDENCE--------------------------------')
mytic=tic;
[Ae,Aee,idxF,idxFx,idxFy,idxFz,Ae1x,Ae1y,Ae1z] = ...
    incidence_matrix3(Kt,[L M N],idxV);
if capacitive_effects_flag==0
    ind=zeros(L*M*N,1);
    ind(ind_c)=val_c;
    ind=ind(idxV);
    ind=setdiff(1:length(idxV),find(ind));
    Aee(max(ind),:)=[];
    idxV(max(ind))=[];
    K = length(idxV); %number of non-empty voxels
end
disp([' Number of DoFs: ', num2str(size(Aee,1)+size(Aee,2))])
disp([' Time for computing incidence ::: ' ,num2str(toc(mytic))]);
disp(' ')
%% Forcing Term: Incident E field
% NOTE: Electric field with components (-iwy/2,iwx/2,0)
% since component x (y) does not depend on x(y), field is calculated at voxel
% barycenter; but in general must be calculated in barycenters of faces!!!
% Thus, we are introducing an approximation here. 
if ext_field_flag
    Ex = Ex_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),omega);
    Ey = Ey_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),omega);
    Ez = Ez_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),omega);
    %%RHS array: <Einc,f_a>_V = int_V dot(Einc,f_a) dV
    Gram = dx*dy*dz; %volume of cubic element
    Vx = (Gram.*Ex)./(dy*dz);
    Vy = (Gram.*Ey)./(dx*dz);
    Vz = (Gram.*Ez)./(dx*dy);
    clear Ex Ey Ez
else
    Vx=zeros(L*M*N,1);
    Vy=zeros(L*M*N,1);
    Vz=zeros(L*M*N,1);
end
%% Matrices Z_real and Z_imag
rho_eF=0.5*(abs(Ae(:,:)).'*rho_eV(:)); clear Ae rho_eV
z_realF=rho_eF;
indFneq=setdiff([1:3*Kt].',idxF);
z_realF(indFneq,:)=0;
z_realx=zeros(L,M,N);
z_realx(idxFx)=z_realF(idxFx);
z_realy=zeros(L,M,N);
z_realy(idxFy)=z_realF(Kt+idxFy);
z_realz=zeros(L,M,N);
z_realz(idxFz)=z_realF(2*Kt+idxFz);
%% Compute Green Tensor 
disp('----COMPUTING GREEN TENSOR--------------------------------')
mytic_G=tic;
[Gmn] = computeGREEN(d,L,M,N,Integration_flag);
disp([' Time for getting Green tensor ::: ' ,num2str(toc(mytic_G))]);
disp(' ')
%% Compute Exponential Tensor For Retardation
if retardation_flag
lambda = co./freq;
ko = 2*pi./lambda;
[Emn] = expGREEN(ko,d,L,M,N);
Gmn = Gmn.*Emn; 
clear Emn %clear exponential tensor
end
%% Compute Circulant Tensors
disp('----COMPUTING CIRCULANT TENSOR--------------------------------')
disp(' Circulant Tensors related to P,L matrices')
mytic_cir=tic;
if capacitive_effects_flag
    [opCirculantP_all,st_sparse_preconP] = computeCIRCULANT(Gmn,d,'P');
else
    opCirculantP_all=[];
end
[opCirculantL_all,st_sparse_preconL] = computeCIRCULANT(Gmn,d,'L');
%%Add constants to Circulants
if capacitive_effects_flag 
    opCirculantP_all = opCirculantP_all/eo;
    st_sparse_preconP = st_sparse_preconP/eo;
else
    st_sparse_preconP=[];
end
opCirculantL_all = (1j*omega*mu)*opCirculantL_all;
st_sparse_preconL = (1j*omega*mu)*st_sparse_preconL;
disp([' Time for getting circulant tensors ::: ' ,num2str(toc(mytic_cir))])
clear Gmn %Green tensor is not used anymore
disp(' ')
%% Generating RHS vector
num_node = size(Aee,1); %all potential nodes in non-empty voxels 
num_curr = size(Aee,2); %all currens in non-empty voxels 
%%Define RHS: (injected currents)
iinj=zeros(L*M*N,1);
iinj(ind_c)=val_c;
iinj=iinj(idxV);
if capacitive_effects_flag
    QIn = zeros(L,M,N);  
    QIn(idxV) = iinj;
    [LfN, MfN, NfN] = size(opCirculantP_all);
    fJ = fftn(QIn(:,:,:),[LfN, MfN, NfN]);
    Jout = opCirculantP_all(:,:,:) .* fJ; 
    JOut = ifftn(Jout);
    JOut = JOut(1:L,1:M,1:N);
    JOut = JOut(idxV);
    if size(JOut,2)>1; JOut=JOut.'; end
    rhs_vect = [Vx(idxFx);Vy(idxFy);Vz(idxFz);-JOut]; 
else
    rhs_vect = [Vx(idxFx);Vy(idxFy);Vz(idxFz);-iinj]; 
end
clear Vx Vy Vz
%% Computing Preconditioner
disp('----COMPUTING PRECONDITIONER--------------------------------')
mytic_prec=tic;
[Y_inv,P_diag,D_diag,LL,UU,PP,QQ,RR] = preparePREC_NEW(d,z_realF,idxFx,idxFy,idxFz,st_sparse_preconP,st_sparse_preconL,Aee,Kt,freq,capacitive_effects_flag);
fPMV = @(JOut_full_in)multiplyPREC_CAP_NEW(JOut_full_in,Aee,Y_inv,P_diag,LL,UU,PP,QQ,RR);
disp([' Time for computing preconditioner ::: ' ,num2str(toc(mytic_prec))]);
disp(' ')
%% Solution of Linear System
disp('----SOLVING LINEAR SYSTEM-------------------------------')
fMVM = @(J) multiplyMATVECT_EDDY_CAP(J,opCirculantL_all,opCirculantP_all,z_realx,z_realy,z_realz,idxF,d,Aee,L,M,N,idxV,freq,capacitive_effects_flag);
mytic_solver=tic;
[vsol] = pgmres_mod(@(J)fMVM(J),rhs_vect, inner_it, tol, outer_it, @(JOut_full_in)fPMV(JOut_full_in) );  
disp([' Time for solving system with gmres ::: ' ,num2str(toc(mytic_solver))]);
disp(' ')
%% extract solution
Jout = zeros(L,M,N,3);
Jout(idxF) = vsol(1:num_curr) ; % return to global variables
%%
%% POST PROCESSING
%%
%% Post Processing J
disp('----POST PROCESSING J------------------------------')
mytic_prec=tic;
[J,XYZ] = fun_my_postRT2(Jout,Kt,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d);
potval=zeros(Kt,1);
potval(idxV)=vsol(num_curr+1:end);
disp([' Total time for post processing J ::: ' ,num2str(toc(mytic_prec))]);
disp(' ')
%% impedance (note that the injected current must be 1A!)
Z=potval(ind_c(1))-potval(ind_c(2)); % Z=V/I
Plosses=0.5*vsol(1:num_curr)'*([z_realx(idxFx)*dx/(dy*dz);z_realy(idxFy)*dy/(dx*dz);z_realz(idxFz)*dz/(dy*dx)].*vsol(1:num_curr)); % losses
Ren=Plosses*2;
disp('-------------------------------------------------------------------')
disp(' ')
if retardation_flag
    disp(['  Impedance: ' ,num2str(Z),' Ohm ']);
else
    disp(['  Impedance: ' ,num2str(Ren+1j*imag(Z)),' Ohm ']);
end
disp(' ')
disp('-------------------------------------------------------------------')
disp(' ')
%% Plot Vectors
if plot_vectorsJ_flag
jjR = real(J);
figure
subplot(1,2,1)
normJR=sqrt(jjR(:,1).^2+jjR(:,2).^2+jjR(:,3).^2);
quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),jjR(:,1),jjR(:,2),jjR(:,3),...
          normJR,4);
axis equal
c1=colorbar;
caxis([min(normJR) max(normJR)]);
xlabel('x')
ylabel('y')
zlabel('z')
title('Current Density Vector \Re Part')
c1.Location = 'southoutside';
xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
%
jjI = imag(J);
subplot(1,2,2)
normJI=sqrt(jjI(:,1).^2+jjI(:,2).^2+jjI(:,3).^2);
quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),jjI(:,1),jjI(:,2),jjI(:,3),...
          normJI,4);
axis equal
c1=colorbar;
caxis([min(normJI) max(normJI)]);
xlabel('x')
ylabel('y')
zlabel('z')
title('Current Density Vector \Im Part')
c1.Location = 'southoutside';
xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
end
%% plot Potential
if plot_potential_flag
xdp = xyz(:,:,:,1);
ydp = xyz(:,:,:,2);
zdp = xyz(:,:,:,3);
figure
subplot(1,2,1)
scatter3(xdp(idxV),ydp(idxV),zdp(idxV),10,'filled','cdata',(real(potval(idxV))))
axis equal
view(3)
colormap jet 
c1=colorbar;
title('Potential \Re')
c1.Location = 'southoutside';
xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
subplot(1,2,2)
scatter3(xdp(idxV),ydp(idxV),zdp(idxV),10,'filled','cdata',imag(potval(idxV)))
axis equal
view(3)
colormap jet 
c2=colorbar;
title('Potential \Im')
c2.Location = 'southoutside';
xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
end
%% paraview
if paraview_export_flag
disp('----EXPORT TO PARAVIEW------------------------------')
xd=xyz(:,:,:,1);
yd=xyz(:,:,:,2);
zd=xyz(:,:,:,3);
xidx=xd(idxV);
yidx=yd(idxV);
zidx=zd(idxV);
P0=[...
    [xidx-dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx+dz/2]];
VP=[1:K;...
    K+1:2*K;...
    2*K+1:3*K;...
    3*K+1:4*K;...
    4*K+1:5*K;...
    5*K+1:6*K;...
    6*K+1:7*K;...
    7*K+1:8*K];
warning off
[~] = ...
    fun_for_ParaView_vec_HEXA(...
    jjR(idxV,:),jjI(idxV,:),P0,VP,dad,[modelname,'J']);
[~] = ...
    fun_for_ParaView_sca_HEXA(...
    real(potval(idxV,:)),imag(potval(idxV,:)),P0,VP,dad,[modelname,'p']);
warning on
end
%%