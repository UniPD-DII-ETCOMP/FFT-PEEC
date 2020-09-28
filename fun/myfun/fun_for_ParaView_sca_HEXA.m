function [Ricc] = fun_for_ParaView_sca_HEXA(sca_r,sca_i,Matrix_P0,VP,dad,name)
% N_volu=size(VP,2);
if exist('res_para','dir') == 0
mkdir('res_para')
end
cd('res_para')
% name='Jc_real';   
disp(['... ',name,' ...'])
[~] = fun_create_xmf_and_hdf_for_ParaView_sca_HEXA(Matrix_P0,VP,sca_r,sca_i,name);
Ricc=1;
cd(dad)
end
function [a1] = fun_create_xmf_and_hdf_for_ParaView_sca_HEXA(Matrix_P0,VP,ssca_re,ssca_im,name)
N_node=size(Matrix_P0,1);
N_hexa=size(VP,2);
delete([name,'.xmf'])
delete([name,'.hdf'])
a1 = floor(N_node/4);
if a1 <=0 
    a1 = 1;
end
a2 = floor(N_hexa/16);
if a2 <=0 
    a2 = 1;
end
aa = name;
a = [aa,'.xmf'];
fileID = fopen(a,'w');
fprintf(fileID,['<?xml version="1.0" ?>','\n']);
fprintf(fileID,['<Xdmf>','\n']);
fprintf(fileID,['<Domain>','\n']);
fprintf(fileID,['<Grid Name="Mesh_Grid">','\n']);
fprintf(fileID,[['<Topology NumberOfElements="',num2str(N_hexa),'" TopologyType="Hexahedron">'],'\n']);
fprintf(fileID,['<DataItem DataType="Int" Dimensions="',num2str(N_hexa),' 8" Format="HDF">','\n']);
fprintf(fileID,[[name,'.hdf:/mesh/hex'],'\n']);
fprintf(fileID,['</DataItem>','\n']);
fprintf(fileID,['</Topology>','\n']);
fprintf(fileID,['<Geometry GeometryType="XYZ">','\n']);
fprintf(fileID,[['<DataItem DataType="Float" Dimensions="',num2str(N_node),' 3" Format="HDF">'],'\n']);
fprintf(fileID,[[name,'.hdf:/mesh/coord'],'\n']);
fprintf(fileID,['</DataItem>','\n']);
fprintf(fileID,['</Geometry>','\n']);
fprintf(fileID,['<Attribute Center="Cell" Name="Scalar_re">','\n']);
fprintf(fileID,[['<DataItem DataType="Float" Dimensions="',num2str(N_hexa),'" Format="HDF">'],'\n']);
fprintf(fileID,[[name,'.hdf:/mesh/sca_re'],'\n']);
fprintf(fileID,['</DataItem>','\n']);
fprintf(fileID,['</Attribute>','\n']);
fprintf(fileID,['<Attribute Center="Cell" Name="Scalar_im">','\n']);
fprintf(fileID,[['<DataItem DataType="Float" Dimensions="',num2str(N_hexa),'" Format="HDF">'],'\n']);
fprintf(fileID,[[name,'.hdf:/mesh/sca_im'],'\n']);
fprintf(fileID,['</DataItem>','\n']);
fprintf(fileID,['</Attribute>','\n']);
fprintf(fileID,['</Grid>','\n']);
fprintf(fileID,['</Domain>','\n']);
fprintf(fileID,['</Xdmf>','\n']);
fclose(fileID);
h5create([name,'.hdf'], '/mesh/coord', size(Matrix_P0.'),'Datatype','double','ChunkSize',[1 a1],'Deflate',4,'FillValue',0.0)
h5write([name,'.hdf'], '/mesh/coord', Matrix_P0.')
h5create([name,'.hdf'], '/mesh/hex', size(VP),'Datatype','int64','ChunkSize',[1 a2],'Deflate',4,'FillValue',int64(0))
h5write([name,'.hdf'], '/mesh/hex', int64(VP-1))
h5create([name,'.hdf'], '/mesh/sca_re', N_hexa,'Datatype','double','ChunkSize',a2,'Deflate',4,'FillValue',0.0)
h5write([name,'.hdf'], '/mesh/sca_re', ssca_re)
h5create([name,'.hdf'], '/mesh/sca_im', N_hexa,'Datatype','double','ChunkSize',a2,'Deflate',4,'FillValue',0.0)
h5write([name,'.hdf'], '/mesh/sca_im', ssca_im)
end