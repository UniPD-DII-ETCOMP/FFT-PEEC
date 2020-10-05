%%
% This function has been inspired by 
% lse_sparse_precon_prepare.m which is available at
% https://github.com/acyucel/VoxHenry and in the directory VoxHenry_functions
%%
function [Y_inv,P_diag,D_diag,LL,UU,PP,QQ,RR] =  preparePREC_NEW(d,...
    z_realF,idxFx,idxFy,idxFz,st_sparse_preconP,...
    st_sparse_preconL,AeeR,Aee,Kt,freq)
omega = 2*pi*freq;
dx = d(1); dy = d(2); dz = d(3);
%%Number of unknowns
num_node = size(Aee,1); % tutti, anche quelli a potenziale vincolato
num_nodeR = size(AeeR,1); % senza quelli a potenziale vincolato
num_curr = size(Aee,2);
num_curr_one3rd = num_curr/3; %since we have 3 basis functions for current
num_currx=length(idxFx);
num_curry=length(idxFy);
num_currz=length(idxFz);
% get the actual values of 'OneoverSigma_e'. This allows different conductivites for each voxel.
diag_pulse=zeros(num_curr,1);
diag_pulse(1:num_currx,1)                                        =1./(z_realF(idxFx)*dx/(dy*dz) + st_sparse_preconL(1));
diag_pulse(num_currx+1:num_currx+num_curry,1)                    =1./(z_realF(Kt+idxFy)*dy/(dz*dx) + st_sparse_preconL(2));
diag_pulse(num_currx+num_curry+1:num_currx+num_curry+num_currz,1)=1./(z_realF(Kt+Kt+idxFz)*dz/(dx*dy) + st_sparse_preconL(3));
%%Sparse 'Y_inv' formation
inds=zeros(num_curr,3);
% rows for sparse 'Y_inv' formation
inds(1:num_curr,1)=[1:1:num_curr];
% columns for sparse 'Y_inv' formation
inds(1:num_curr,2)=inds(1:num_curr,1);
% values for sparse 'Y_inv' formation, split in three sets:
% first set is from 1 to 3/5 of num_curr, i.e. Ix, Iy, Iz
inds(:,3)=(diag_pulse);
%%Create Sparse 'Y_inv'
Y_inv=sparse(inds(:,1),inds(:,2),inds(:,3));
%%Obtain Schur complement   
% [A  B] = [Z       Ae']
% [C  D]   [Ae        0]
Sch_comp= D_diag - P_diag*(AeeR*Y_inv*AeeR.');        
%%Computing Schur Inversion
[LL,UU,PP,QQ,RR] = lu(Sch_comp);               
end
