% % r = table2array(re);
% result = string(zeros(length(r),1));
% for i = 1:length(r)
%     result(i) = strcat('E',r(i));
% end
% filename = 're.xlsx';
% xlswrite(filename,result,'sheet1','A');  %��EXCEL�ļ��� T �н��ж�ȡ
excelname = 'E:\Design\Carotid_Analysis\Bifurcation_Task\����\re.xlsx';
[~,filename,~ ]= xlsread(excelname,'sheet1','A:A');  %��EXCEL�ļ��� A �н��ж�ȡ
filename = filename(2:end);
input_bvdat = string(zeros(length(filename),1));
input_centerlinedat = string(zeros(length(filename),1));
inputdir = 'E:\Design\Carotid_Analysis\Bifurcation_Task\����\data'; % ����
for i = 1:length(filename)
    %     input_bvdat{i} = strcat(filename{i},'_bv.dat');
    %     input_centerlinedat{i} = strcat(filename{i},'_centerlines.dat');
    input_bvdat(i) = strcat(inputdir,'\',filename{i},'_bv.dat');
    input_centerlinedat(i) = strcat(inputdir,'\',filename{i},'_centerlines.dat');
end
%%%�ļ������ļ����У�exist����ֵΪ2������Ϊ0��
for i = 1:length(filename)
    filebv(i) = exist(input_bvdat(i),'file');
    filecen(i) = exist(input_centerlinedat(i),'file');
end
%�� ��
[m,n]=find(filebv==0)
[m,n]=find(filecen==0)
