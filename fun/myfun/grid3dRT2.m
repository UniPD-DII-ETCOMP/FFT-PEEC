function r = grid3dRT2(x, y, z)
L = length(x);
M = length(y);
N = length(z);
r = zeros(L,M,N,3);
try 
parfor ix = 1:L
    xx = x(ix);
    for iy = 1:M
        yy = y(iy);
        for iz = 1:N
            zz = z(iz);
            r(ix,iy,iz,:) = [xx yy zz];
        end
    end
end
catch
for ix = 1:L
    xx = x(ix);
    for iy = 1:M
        yy = y(iy);
        for iz = 1:N
            zz = z(iz);
            r(ix,iy,iz,:) = [xx yy zz];
        end
    end
end    
end
end
