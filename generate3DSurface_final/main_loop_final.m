d = dir( angleDir );
for i = 1:length( d )
   if ( strcmp( d( i ).name, '.' ) || strcmp( d( i ).name, '..' ) || strcmp( d( i ).name, 'E38162_R_20110127_094433' ) )
       continue;
   end

   if ( ~d( i ).isdir )
       continue;
   end
   
   disp( [ 'on ', int2str( i ), ' of ', int2str( length( d ) ) ] );

   idx = strfind( d( i ).name, '_' );
   proj = d( i ).name( 1:idx(2)-1 );
   fullDir = [ angleDir, '/', d( i ).name, '/' ];
   
   projFile = [ proj, '.QVJ' ];
   
   fullDir = strrep( fullDir, '\', '/' );
   tmpDir = strrep( fullDir, '/', '\' );
   [centers, newcutval]=Cascade3DOutput( tmpDir, projFile, [ outDir, '/', proj ], sm, false, majorweighting,islumen,cutoff);
  
   cx = num2str( centers( 1, 1 ) );
   cy = num2str( centers( 1, 2 ) );
   cz = num2str( centers( 1, 3 ) );
   
   ix = num2str( centers( 2, 1 ) );
   iy = num2str( centers( 2, 2 ) );
   iz = num2str( centers( 2, 3 ) );
   
   ex = num2str( centers( 3, 1 ) );
   ey = num2str( centers( 3, 2 ) );
   ez = num2str( centers( 3, 3 ) );
   
   geoCmd = [ 'vmtksurfacereader -ifile "', outDir, '/', proj, '.vtk" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkcenterlines -seedselector "pointlist" -sourcepoints ', cx, ' ', cy, ' ', cz, ' -targetpoints ', ix, ' ', iy, ' ', iz, ' ', ex, ' ', ey, ' ', ez, ' -ofile "', outDir, '/', proj, '.vtp" ' ]; 
   geoCmdArray{ end + 1 } = geoCmd;
   geoCmd = [ 'vmtkcenterlineattributes -ifile "', outDir, '/', proj, '.vtp" --pipe vmtkbranchextractor -radiusarray@ MaximumInscribedSphereRadius --pipe vmtkbifurcationreferencesystems --pipe vmtkbifurcationvectors -ofile "', outDir, '/', proj, '_bv.dat" ' ];
   geoCmdArray{ end + 1 } = geoCmd;
   geoCmd = [ 'vmtkcenterlineattributes -ifile "', outDir, '/', proj, '.vtp" --pipe vmtkbranchextractor -radiusarray@ MaximumInscribedSphereRadius --pipe vmtkbranchgeometry -ofile "', outDir, '/', proj, '_clcg.dat" ' ];
   geoCmdArray{ end + 1 } = geoCmd;
      
   viewCmd = [ 'vmtksurfacereader -ifile "', outDir, '/', proj, '.vtk" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtkrenderer --pipe vmtksurfaceviewer -opacity 0.25 --pipe vmtksurfaceviewer -ifile "', outDir, '/', proj, '.vtp" -array MaximumInscribedSphereRadius' ];
   viewCmdArray{ end + 1 } = viewCmd;
   
   viewCmd2 = [ 'vmtksurfacereader -ifile "', outDir, '/', proj, '.vtk" --pipe vmtksurfacesmoothing -passband 0.1 -iterations 30 --pipe vmtksurfaceviewer -opacity 0.25' ];
   viewCmdArray2{ end + 1 } = viewCmd2;
   
   areaCmd = [ 'START /WAIT CMD /C area_calc2.bat "', outDir, '/', proj, '"' ];
   areaCmdArray{ end + 1 } = areaCmd;
   
   tc24 = [ 'START /WAIT CMD /C tort_calc_24.bat "', outDir, '/', proj, '"' ];
   tc25 = [ 'START /WAIT CMD /C tort_calc_25.bat "', outDir, '/', proj, '"' ];
   tc35 = [ 'START /WAIT CMD /C tort_calc_35.bat "', outDir, '/', proj, '"' ];
   tortCmdArray24{ end + 1 } = tc24;
   tortCmdArray25{ end + 1 } = tc25;
   tortCmdArray35{ end + 1 } = tc35;
   
   cutoffArray{ end + 1 } = newcutval;
end

disp('done!');

