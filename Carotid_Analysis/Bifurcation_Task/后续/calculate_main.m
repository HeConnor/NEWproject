function calculate_main(excelname,inputdir)
% excelname = 'E:\Design\Carotid_Analysis\Bifurcation_Task\����\case1_ruolan_sm_100_����.xlsx';
% outputexcel = 'E:\Design\Carotid_Analysis\Bifurcation_Task\����\resu279.xlsx';
[~,filename,~ ]= xlsread(excelname,'sheet1','A:A');  %��EXCEL�ļ��� A �н��ж�ȡ
filename = filename(2:end);
input_bvdat = string(zeros(length(filename),1));
input_centerlinedat = string(zeros(length(filename),1));
%%%����Ҫ�ж���������Դ�Ƿ����ƶ��ļ����У�����
% inputdir = 'E:\Design\Carotid_Analysis\case\case1\ruolan_outdir\sm_100'; %%%����
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
