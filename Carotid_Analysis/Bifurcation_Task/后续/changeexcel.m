% % r = table2array(re);
% result = string(zeros(length(r),1));
% for i = 1:length(r)
%     result(i) = strcat('E',r(i));
% end
% filename = 're.xlsx';
% xlswrite(filename,result,'sheet1','A');  %对EXCEL文件第 T 列进行读取
excelname = 'E:\Design\Carotid_Analysis\Bifurcation_Task\后续\re.xlsx';
[~,filename,~ ]= xlsread(excelname,'sheet1','A:A');  %对EXCEL文件第 A 列进行读取
filename = filename(2:end);
input_bvdat = string(zeros(length(filename),1));
input_centerlinedat = string(zeros(length(filename),1));
inputdir = 'E:\Design\Carotid_Analysis\Bifurcation_Task\后续\data'; % 更改
for i = 1:length(filename)
    %     input_bvdat{i} = strcat(filename{i},'_bv.dat');
    %     input_centerlinedat{i} = strcat(filename{i},'_centerlines.dat');
    input_bvdat(i) = strcat(inputdir,'\',filename{i},'_bv.dat');
    input_centerlinedat(i) = strcat(inputdir,'\',filename{i},'_centerlines.dat');
end
%%%文件存在文件夹中，exist返回值为2，否则为0；
for i = 1:length(filename)
    filebv(i) = exist(input_bvdat(i),'file');
    filecen(i) = exist(input_centerlinedat(i),'file');
end
%行 列
[m,n]=find(filebv==0)
[m,n]=find(filecen==0)
