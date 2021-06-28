function [ col ] = hexa_mesh2(VP,P0,alpha,col)

patch('Faces',VP(1:4,:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha);%,'edgecolor',[0.7 0.7 0.7])
hold on
patch('Faces',VP(5:8,:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha);%'edgecolor',[0.7 0.7 0.7])
patch('Faces',VP([1,2,6,5],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha);%,'edgecolor',[0.7 0.7 0.7])
patch('Faces',VP([4,3,7,8],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha);%,'edgecolor',[0.7 0.7 0.7])
patch('Faces',VP([1,4,8,5],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha);%,'edgecolor',[0.7 0.7 0.7])
patch('Faces',VP([2,3,7,6],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha);%,'edgecolor',[0.7 0.7 0.7])


end

