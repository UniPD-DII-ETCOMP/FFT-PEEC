%%
% This function has been inspired by (a part of)
% lse_generate_circulant_tensor.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [Gp_mn] = circulant(G_mn,NUM)
[L, M, N, ~] = size(G_mn);
Gp_mn = zeros(2*L,2*M,2*N,NUM);
% Cube 'L'
[Gp_coeff_L] = gperiodic_coeff_nop('L');
% Cube 'M'
[Gp_coeff_M] = gperiodic_coeff_nop('M');
% Cube 'N'
[Gp_coeff_N] = gperiodic_coeff_nop('N');
% Cube 'LM'
[Gp_coeff_LM] = gperiodic_coeff_nop('LM');
% Cube 'LN'
[Gp_coeff_LN] = gperiodic_coeff_nop('LN');
% Cube 'MN'
[Gp_coeff_MN] = gperiodic_coeff_nop('MN');
% Cube 'LMN' equal to 1, no calculate [Gp_coeff_LMN] = gperiodic_coeff('LMN');
Gp_mn(1:L,1:M,1:N,:) = G_mn; 
for ii = 1:NUM
    % Cube 'L'
    Gp_mn(L+2:2*L,1:M,1:N,ii)         = G_mn(L:-1:2,1:M,1:N,ii) * Gp_coeff_L(ii,1);
    % Cube 'M'
    Gp_mn(1:L,M+2:2*M,1:N,ii)         = G_mn(1:L,M:-1:2,1:N,ii) * Gp_coeff_M(ii,1);
    % Cube 'N'
    Gp_mn(1:L,1:M,N+2:2*N,ii)         = G_mn(1:L,1:M,N:-1:2,ii) * Gp_coeff_N(ii,1);
    % Cube 'LM'
    Gp_mn(L+2:2*L,M+2:2*M,1:N,ii)     = G_mn(L:-1:2,M:-1:2,1:N,ii) * Gp_coeff_LM(ii,1);
    % Cube 'LN'
    Gp_mn(L+2:2*L,1:M,N+2:2*N,ii)     = G_mn(L:-1:2,1:M,N:-1:2,ii) * Gp_coeff_LN(ii,1);
    % Cube 'MN'
    Gp_mn(1:L,M+2:2*M,N+2:2*N,ii)     = G_mn(1:L,M:-1:2,N:-1:2,ii) * Gp_coeff_MN(ii,1);
end
% Cube 'LMN'
Gp_mn(L+2:2*L,M+2:2*M,N+2:2*N,:) = G_mn(L:-1:2,M:-1:2,N:-1:2,:) ;
end