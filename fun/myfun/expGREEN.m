function [E] = expGREEN(k,d,L,M,N)
E = zeros(L,M,N,1);
% infofN = whos('E'); memestimated = 2*infofN.bytes/(1024*1024);
% disp(['  Memory for temporarily storing Exp Tensor (MB) ::: ' , num2str(memestimated)]);
% disp('  Start computing E tensor')
tic
for mx = 1:L
    for my = 1:M
        for mz = 1:N
                m = [mx,my,mz];
                xyz_trg = ((m-1) .* d)';
                E(mx,my,mz,:) = exp(-1i*k*norm(-xyz_trg));
        end
    end
end
clear xyz_src xyz_trg   
% disp(['  Time for exponetial factor green tensor ::: ',num2str(toc)])
end

