%%
% This function has been inspired by (a part of)
% lse_generate_circulant_tensor.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [G] = computeGREEN(d,L,M,N,Integration)
NGx = 2; %number of integration points x-dir
NGy = 2; %number of integration points y-dir
NGz = 2; %number of integration points z-dir
G = zeros(L,M,N,1); %1 because we have 3 bases function but without considering normalization, the tensors are the same
if strcmp(Integration,'NumAn') 
[face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross] = prepareINTEGRATION(d);
[xyz_src,ww] = GaussPoints([NGx,NGy,NGz],d,Integration);
try 
   G=computeGREEN_f90_mexed(L,M,N,xyz_src,length(ww),face,r_f,r_e,u_e,triang_face,idxf_dot,...
                            idxf_cross,ww,d.',22); 
    G=reshape(G,L,M,N);
catch
    warning('============================================================')
    warning('- MEX function not supported, try to re-mex it: run /fortran/make.m.')
    disp(' ')
    warning('- Switched to slow Matlab function.') 
    disp(' ')
    warning('- You can use NumNum instead (much faster, less accuracy)')
    warning('============================================================')
    for mx = 1:L
            for my = 1:M
                for mz = 1:N
                    m = [mx,my,mz];
                    r_m = (m-1) .* d;
                    xyz_trg = xyz_src + r_m; %Gauss points on target voxel
                    G(mx,my,mz,:) =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,ww);
                end
            end
    end
end
clear face face_diff r_f r_e u_e triang_face idxf_dot idxf_cross
clear xyz_src xyz_trg w
elseif strcmp(Integration,'NumNum')
[xyz1,ww1,xyz2,ww2] = GaussPoints([NGx,NGy,NGz],d,Integration);
Lii = 0;
for ii = 1:NGx*NGy*NGz
    for jj = 1:(NGx+1)*(NGy+1)*(NGz+1)
       Lii = Lii+ww1(ii)*ww2(jj)/(sqrt(sum((xyz1(ii,:)-xyz2(jj,:)).^2)));
    end
end
Lii = Lii/(4*pi);
    parfor mx = 1:L
        for my = 1:M
            for mz = 1:N
                m = [mx,my,mz];
                xyz_trg = ((m-1) .* d)';
                G(mx,my,mz,:) = Integrate_NumNum(xyz_trg,d,Lii);
            end
        end
    end
    
clear xyz1 ww1 xyz2 ww2
clear xyz_trg
end   
end

