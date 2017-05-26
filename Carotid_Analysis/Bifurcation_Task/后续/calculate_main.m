function calculate_main(excelname,inputdir)
% excelname = 'E:\Design\Carotid_Analysis\Bifurcation_Task\后续\case1_ruolan_sm_100_副本.xlsx';
% outputexcel = 'E:\Design\Carotid_Analysis\Bifurcation_Task\后续\resu279.xlsx';
[~,filename,~ ]= xlsread(excelname,'sheet1','A:A');  %对EXCEL文件第 A 列进行读取
filename = filename(2:end);
input_bvdat = string(zeros(length(filename),1));
input_centerlinedat = string(zeros(length(filename),1));
%%%首先要判断所需数据源是否在制定文件夹中！！！
% inputdir = 'E:\Design\Carotid_Analysis\case\case1\ruolan_outdir\sm_100'; %%%更改
for i = 1:length(filename)
    %     input_bvdat{i} = strcat(filename{i},'_bv.dat');
    %     input_centerlinedat{i} = strcat(filename{i},'_centerlines.dat');
    input_bvdat(i) = strcat(inputdir,'\',filename{i},'_bv.dat');
    input_centerlinedat(i) = strcat(inputdir,'\',filename{i},'_centerlines.dat');
end
zout = {};
for i = 1:length(input_bvdat)
    [results,pointradius] = findpoint(input_bvdat(i),input_centerlinedat(i));
    if i==1
        zout = [results,pointradius];
    else
        zout = [zout;results(2,:),pointradius(2,:)];
    end
end

% zout = [['filename';filename],zout];
xlswrite(excelname,zout,'sheet1','B');
disp('done!');
end
