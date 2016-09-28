clear all;

%%%% Change to path with exams downloaded from VTS
% angleDir = 'C:/Users/Dan Admin/Desktop/2010-10-06_angles/DHS';
angleDir = 'E:/Design/E99078_L_20160421_182139';
% try 3 different amounts of smoothing
% make separate outDir for each level of smoothing
% run one at a time, using comment symbols (%) to select them
 
sm = 100;
% outDir = 'C:/Users/Dan Admin/Desktop/2010-10-06_angles/DHS_OUT_100';
outDir = 'E:/Design/generate3DSurface_final/OUT_dir';

%sm = 200;
%outDir = 'C:/Users/Dan Admin/Desktop/2010-10-06_angles/DHS_OUT_200';

%sm = 300;
%outDir = 'C:/Users/Dan Admin/Desktop/2010-10-06_angles/DHS_OUT_300';


% do not change
majorweighting = 2;
islumen = true;
cutoff = -1;

geoCmdArray = {};
viewCmdArray = {}; % for viewing surface + centerlines
viewCmdArray2 = {}; % just for viewing surface
areaCmdArray = {};
tortCmdArray24 = {};
tortCmdArray25 = {};
tortCmdArray35 = {};
cutoffArray = {};

main_loop_final

fid = fopen( [ outDir, '/geo_cmds.txt' ], 'w' );
for i=1:length( geoCmdArray )
    fprintf(fid, '%s \n\n', geoCmdArray{ i } );
end
fclose( fid );

fid = fopen( [ outDir, '/view_cmds.txt' ], 'w' );
for i=1:length( viewCmdArray )
    fprintf(fid, '%s \n\n', viewCmdArray{ i } );
end
fclose( fid );

fid = fopen( [ outDir, '/view_surface_cmds.txt' ], 'w' );
for i=1:length( viewCmdArray2 )
    fprintf(fid, '%s \n\n', viewCmdArray2{ i } );
end
fclose( fid );

fid = fopen( [ outDir, '/area_cmds.txt' ], 'w' );
for i=1:length( areaCmdArray )
    fprintf(fid, '%s \n\n', areaCmdArray{ i } );
end
fclose( fid );

fid = fopen( [ outDir, '/tort_cmds.txt' ], 'w' );
for i=1:length( tortCmdArray24 )
    fprintf(fid, '%s \n\n', tortCmdArray24{ i } );
    fprintf(fid, '%s \n\n', tortCmdArray25{ i } );
    fprintf(fid, '%s \n\n', tortCmdArray35{ i } );
end
fclose( fid );