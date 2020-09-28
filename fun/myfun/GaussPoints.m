function [xyz1,ww1,varargout] = GaussPoints(N,d,type)
%%Warmup Procedure
dx = d(1); dy = d(2); dz = d(3);
%%Since a set of points is required also in numerical/numerical, calcolate
%%one directly.
[x,wx]=lgwt(N(1),-dx/2,dx/2);
[y,wy]=lgwt(N(2),-dy/2,dy/2);
[z,wz]=lgwt(N(3),-dz/2,dz/2);
%From line to cube 
xyz1=zeros(N(1)*N(2)*N(3),3);
ww1=zeros(N(1)*N(2)*N(3),1);
hh=1;
for ii = 1:N(1)
    for jj = 1:N(2)
        for kk = 1:N(3)
            xyz1(hh,:)=[x(ii) y(jj) z(kk)];
            ww1(hh)=wx(ii)*wy(jj)*wz(kk);
            hh=hh+1;
        end
    end
end
if strcmp(type,'NumNum') %numerical/numerical
    %Compute a second set of points
    N2 = N+1;
    %%Edge lengths
    [x,wx]=lgwt(N2(1),-dx/2,dx/2);
    [y,wy]=lgwt(N2(2),-dy/2,dy/2);
    [z,wz]=lgwt(N2(3),-dz/2,dz/2);
    %From line to cube 
    xyz2=zeros(N2(1)*N2(2)*N2(3),3);
    ww2=zeros(N2(1)*N2(2)*N2(3),1);
    hh=1;
    for ii = 1:N2(1)
        for jj = 1:N2(2)
            for kk = 1:N2(3)
                xyz2(hh,:)=[x(ii) y(jj) z(kk)];
                ww2(hh)=wx(ii)*wy(jj)*wz(kk);
                hh=hh+1;
            end
        end
    end
    varargout{1} = xyz2;
    varargout{2} = ww2;
end 
end

