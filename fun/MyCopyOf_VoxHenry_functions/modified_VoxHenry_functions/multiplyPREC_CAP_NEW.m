%%
% This function has been inspired by 
% lse_sparse_precon_multiply.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [JOut_full_out] = multiplyPREC_CAP_NEW(JOut_full_in,AeR,Y_inv,P_diag,LL,UU,PP,QQ,RR)
num_nodeR=size(AeR,1); 
num_curr=size(AeR,2);
JOut_full_out=zeros(num_nodeR+num_curr,1);
warning off
JOut_full_out(num_curr+1:num_curr+num_nodeR) = ...
    QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_nodeR) - ( P_diag * (AeR * (Y_inv * JOut_full_in(1:num_curr))) ) )))));
warning on
JOut_full_out(1:num_curr) = Y_inv*(JOut_full_in(1:num_curr) - ((AeR.')*JOut_full_out(num_curr+1:num_curr+num_nodeR)));
end