VP=zeros(8,length(idx));
Lp=L+1;
Mp=M+1;
Np=N+1;
for ii = 1:length(idx)
    lmn=idx2triplet(idx(ii),L,M,N);
    l=lmn(1);
    m=lmn(2);
    n=lmn(3);
    VP(1,ii)=l  +Lp*(m  -1)+Lp*Mp*(n+1-1);
    VP(2,ii)=l  +Lp*(m  -1)+Lp*Mp*(n-1  );
    VP(3,ii)=l+1+Lp*(m  -1)+Lp*Mp*(n-1  );
    VP(4,ii)=l+1+Lp*(m  -1)+Lp*Mp*(n+1-1);
    VP(5,ii)=l  +Lp*(m+1-1)+Lp*Mp*(n+1-1);
    VP(6,ii)=l  +Lp*(m+1-1)+Lp*Mp*(n-1  );
    VP(7,ii)=l+1+Lp*(m+1-1)+Lp*Mp*(n-1  );
    VP(8,ii)=l+1+Lp*(m+1-1)+Lp*Mp*(n+1-1);
end
P0=zeros(max(VP(:)),3);
for ii = 1:max(VP(:))
P0(ii,:)=idx2triplet(ii,L+1,M+1,N+1).*[dx dy dz];
end
figure
hexa_mesh2(VP,P0,0.3,'r')
view(3)
title('mesh')
axis off
axis equal
set(gcf,'color','w');