%%
% This function has been inspired by (a part of)
% lse_generate_circulant_tensor.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [Circ,Prec_Circ] = computeCIRCULANT(G,d,matrix)
[L,M,N] = size(G);
dx = d(1); dy = d(2); dz = d(3);
tic 
[Gp] = circulant(G,1);
Gp_mn = zeros(2*L,2*M,2*N,3);
if strcmp(matrix ,'L')
    Gp_mn(:,:,:,1) = Gp/(dy^2*dz^2); 
    Gp_mn(:,:,:,2) = Gp/(dx^2*dz^2);
    Gp_mn(:,:,:,3) = Gp/(dx^2*dy^2); 
elseif strcmp(matrix ,'P')
    vol = dx*dy*dz;
    volvol = vol*vol;
    Gp_mn = Gp/volvol;
end
if strcmp(matrix ,'L')
    Prec_Circ(1)=squeeze(G(1,1,1,1)/(dy^2*dz^2));
    Prec_Circ(2)=squeeze(G(1,1,1,1)/(dx^2*dz^2));
    Prec_Circ(3)=squeeze(G(1,1,1,1)/(dx^2*dy^2));
elseif strcmp(matrix ,'P')
    Prec_Circ=squeeze(G(1,1,1,1)/volvol);
end
    tic
    Circ = fft_operator(Gp_mn);
end

