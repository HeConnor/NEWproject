excelDir = 'E:\Design\Carotid_Analysis\Bifurcation_Task\ºóÐø\results';
partpath0 = 'E:\Design\Carotid_Analysis\case';
d = dir( excelDir );
for i = 1:length( d )
   if ( strcmp( d( i ).name, '.' ) || strcmp( d( i ).name, '..' )) 
       continue;
   end
   disp( [ 'on ', int2str( i-2 ),'¡ª¡ª',d( i ).name ] );
   idx = strfind( d( i ).name, '_' );
   partpath1 = d( i ).name(1:idx(1)-1); %case x
   partpath2 = strcat(d( i ).name(idx(1)+1:idx(2)),'outdir');
   partpath3 = d( i ).name(idx(2)+1:idx(3)+3); %sm_x00
   dat_filepath = strcat(partpath0,'\',partpath1,'\',partpath2,'\',partpath3);
   excelname = strcat(excelDir,'\',d( i ).name);
   disp([{dat_filepath};{excelname}]);
   calculate_main(excelname,dat_filepath);
end
disp('All done!');

% excelname
% input_bvdat
% input_centerlinedat
