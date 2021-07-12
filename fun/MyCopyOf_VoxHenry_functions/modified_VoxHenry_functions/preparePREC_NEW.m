%%
% This function has been inspired by 
% lse_sparse_precon_prepare.m which is available at
% https://github.com/acyucel/VoxHenry and in the directory VoxHenry_functions
%%
function [Y_inv,P_diag,D_diag,LL,UU,PP,QQ,RR] =  preparePREC_NEW(d,...
    z_realF,idxFx,idxFy,idxFz,st_sparse_preconP,...
    st_sparse_preconL,Aee,Kt,freq,capacitive_effects_flag)
omega = 2*pi*freq;
dx = d(1); dy = d(2); dz = d(3);
num_node = size(Aee,1); 
num_curr = size(Aee,2);
num_currx=length(idxFx);
num_curry=length(idxFy);
num_currz=length(idxFz);
diag_pulse=zeros(num_curr,1);
diag_pulse(1:num_currx,1)                                        =1./(z_realF(idxFx)*dx/(dy*dz) + st_sparse_preconL(1));
diag_pulse(num_currx+1:num_currx+num_curry,1)                    =1./(z_realF(Kt+idxFy)*dy/(dz*dx) + st_sparse_preconL(2));
diag_pulse(num_currx+num_curry+1:num_currx+num_curry+num_currz,1)=1./(z_realF(Kt+Kt+idxFz)*dz/(dx*dy) + st_sparse_preconL(3));
inds=zeros(num_curr,3);
inds(1:num_curr,1)=[1:1:num_curr];
inds(1:num_curr,2)=inds(1:num_curr,1);
inds(:,3)=(diag_pulse);
Y_inv=sparse(inds(:,1),inds(:,2),inds(:,3));
if capacitive_effects_flag
    diag_pulse = zeros(num_node,1); 
    diag_pulse(1:num_node,1) = st_sparse_preconP(1); 
    inds = zeros(num_node,3);
    inds(1:num_node,1) = [1:1:num_node];
    inds(1:num_node,2)=inds(1:num_node,1);
    inds(:,3)=(diag_pulse);
    P_diag =sparse(inds(:,1),inds(:,2),inds(:,3));
    diag_pulse = zeros(num_node,1); 
    diag_pulse(1:num_node,1) = -1i*omega; 
    inds = zeros(num_node,3);
    inds(1:num_node,1) = [1:1:num_node];
    inds(1:num_node,2)=inds(1:num_node,1);
    inds(:,3)= diag_pulse;
    D_diag =sparse(inds(:,1),inds(:,2),inds(:,3));
    % [A  B] = [Z       Ae']
    % [C  D]   [P*Ae  -iwI ]
    Sch_comp= D_diag - P_diag*(Aee*Y_inv*Aee.');
    [LL,UU,PP,QQ,RR] = lu(Sch_comp);
else
    % [A  B] = [Z    Ae']
    % [C  D]   [Ae    0 ]
    P_diag=1;
    D_diag=0;
    Sch_comp=  - (Aee*Y_inv*Aee.');
    [LL,UU,PP,QQ,RR] = lu(Sch_comp);    
end
end
