function [results,pointradius]=findpoint(input_bvdat,input_centerlinedat)
% clear all;
% close all;
% �����excel��һ���ļ�·�����ݣ��Ӷ����Դ�bv��centerlines_dat�ļ�
% b=importdata('D:\Program Files\Matlab\�༭\task1\data\E720_L_bv.dat');
b=importdata(input_bvdat);
x=b.data(1:3,1:3);
%�������������ݣ���ʼ�ҵ�
% a=importdata('D:\Program Files\Matlab\�༭\task1\data\E720_L_centerlines.dat');		%���ļ�
a=importdata(input_centerlinedat);
data=a.data;
% plot3(data(:,1),data(:,2),data(:,3));		%����������ͼ
% hold on;
% plot3(x(:,1),x(:,2),x(:,3),'r+');			%bv.dat����������
N=length(data);   %�������ܳ���
%ICA(1/2/3/4/5)����ȡ+ICA0
C=repmat(x(3,1:3),N,1);			%repmat ���ƺ�ƽ�̾���
%bv.dat��groupID 3 �����������γ���data������ͬ�ľ���
dis=sum((C-data(:,1:3)).^2,2);		%���룬������ÿ�㵽bv����groupID 3���ľ����ƽ����
index=zeros(6,1);
r=index;
[dismin,index(1)]=min(dis);     %ica0 ��ֵ��С�ĵ��λ�ú����ֵ
indexs=find(data(:,8)==0);  %�����֧���ն˱��
stapoint = indexs(1)+1; %�ҳ�ICA֧��ʼ��,��ECAĩ�˿�ʼ��
dis=dis(stapoint:index(1)-1);  %%��ICA1~5����ICA0������ʼ��Ѱ��
%%%disΪǰ100�㵽ICA0�ļ��ƽ��
r(1)=data(index(1),4);      %ICA0�뾶
for i=2:6
      %%%���dis�������100��뾶��ƽ�� ��ֵ
    [dismin,index1]=min(abs(dis-data(stapoint:index(i-1)-1,4).^2));
    if isempty(index1)
        i;r(i:6)=0 ; %�㲻����
        break;
    end
    index(i)=stapoint+index1-1; %%%ת��ΪICA���λ��
    r(i)=data(index(i),4); %%%����뾶
    if (dismin>0.5*r(i))
        i;r(i:6)=0 ;       %�����ж�
        break;
    end
    xx=data(index(i),1:3);  
    NICA=length(data((stapoint:index(i)-1),4));
    C=repmat(xx,NICA,1);  %�ҵ�����һ�㲢���������ǰ50��ľ��룬Ϊ��һ����׼��
%     plot3(xx(1),xx(2),xx(3),'b*');
    dis=sum((C-data((stapoint:index(i)-1),1:3)).^2,2);
end
% [dismin,index(1)]=min(dis);     %ica0
% dis=dis(index(1)-100:index(1)-1);
% r(1)=data(index(1),4);
% for i=2:6
%     [dismin,index1]=min(abs(dis-data(index(i-1)-100:index(i-1)-1,4).^2));
%     index(i)=index(i-1)+index1-101;
%     r(i)=data(index(i),4);
%     if (dismin>0.5*r(i))
%         i;r(i:6)=0 ;       %�����ж�
%         break;
%     end
%     xx=data(index(i),1:3);
%     C=repmat(xx,100,1);
%     %     plot3(xx(1),xx(2),xx(3),'b*');
%     dis=sum((C-data((index(i)-100:index(i)-1),1:3)).^2,2);
% end
ICA=[r index];

%ECA(0/1)����ȡ
index=zeros(2,1);
r=index;
%bv.dat��groupID 2 �����������γ���data������ͬ�ľ���
C=repmat(x(2,1:3),N,1);  %ECA0
dis=sum((C-data(:,1:3)).^2,2);%���룬������ÿ�㵽bv����groupID 2���ľ����ƽ����
[dismin,index(1)]=min(dis);
indexs=index(1)-100; %%����ECA1֮ǰ�ĵ��һ�����
if indexs<=0
    indexs=1;
end
dis=dis(indexs:index(1)-1);
r(1)=data(index(1),4);
for i=2:2
    [dismin,index1]=min(abs(dis-data(indexs:index(i-1)-1,4).^2));
    %index(i)=index(i-1)+index1-201;
    index(i)=indexs+index1-1;
    r(i)=data(index(i),4);
    if (dismin>0.5*r(i))
        i;r(i:2)=0 ;
        break;
    end
    xx=data(index(i),1:3);
    C=repmat(xx,100,1);
%     plot3(xx(1),xx(2),xx(3),'r*');
    if index(i)-100<1
        indexs=1;
        C=repmat(xx,index(i)-1,1);
    else
        indexs=index(i)-100;
    end
    dis=sum((C-data((indexs:index(i)-1),1:3)).^2,2);
end
ECA=[r index];

index=zeros(6,1);
r=index;
indexs=find(data(:,8)==0);  %�����֧���ն˱��
%     ica0=data(ICA(1,2),:);
%     ecao=data(ECA(1,2),:);
%�ڵ�һ֧ECA���ҵ㡪�� ECA -1 -2 -3 -4 -5
N=indexs(1)-ECA(1,2); %ECA0���ն˵���
C=repmat(x(2,1:3),N,1);			%repmat ���ƺ�ƽ�̾���
index(1)=ECA(1,2);
r(1)=data(ECA(1,2),4);
dis=sqrt(sum((C-data(ECA(1,2)+1:indexs(1),1:3)).^2,2));%ECA0-��֧��ECA0����
for i=2:6
    [dismin,index1]=min(abs(dis-data(index(i-1)+1:indexs(1),4)));
    if isempty(index1)
        i;r(i:4)=0 ;
        break;
    end
    index(i)=index(i-1)+index1;
    r(i)=data(index(i),4);
    if (dismin>0.5*r(i))
        i;r(i:6)=0 ;       %����
        break;
    end
    xx=data(index(i),1:3);
    N=length(data((index(i)+1:indexs(1)),4));
    C=repmat(xx,N,1);
%     plot3(xx(1),xx(2),xx(3),'b*');
    dis=sqrt(sum((C-data((index(i)+1:indexs(1)),1:3)).^2,2));
end
ECA_down=[r index];
%�ڵڶ�֧ICA���ҵ㡪�� ICA -1 -2 -3 -4 -5
index=zeros(6,1);
r=index;
N=indexs(2)-ICA(1,2); %ICA0����֧ĩ��
C=repmat(x(3,1:3),N,1);			%repmat ���ƺ�ƽ�̾���
index(1)=ICA(1,2);
r(1)=data(ICA(1,2),4);
dis=sqrt(sum((C-data(ICA(1,2)+1:indexs(2),1:3)).^2,2));
for i=2:6
    [dismin,index1]=min(abs(dis-data(index(i-1)+1:indexs(2),4)));
    if isempty(index1)
        i;r(i:4)=0 ;
        break;
    end
    index(i)=index(i-1)+index1;
    r(i)=data(index(i),4);
    if (dismin>0.5*r(i))
        i;r(i:6)=0 ;       %����
        break;
    end
    xx=data(index(i),1:3);
    N=length(data((index(i)+1:indexs(2)),4));
    C=repmat(xx,N,1);
%     plot3(xx(1),xx(2),xx(3),'b*');
    dis=sqrt(sum((C-data((index(i)+1:indexs(2)),1:3)).^2,2));
end
ICA_down=[r index];
%�ݶ�ΪICA֧�õ�CCA�����е�
%��������ICA_down��ECA_down�ж�Ӧ����һ��Ϊ0����֧ͬʱ�������ݣ����߶�������
index=zeros(4,1);
r=index;
index(1)=ICA_down(3,2);
r(1)=ICA_down(3,1);     %�ٶ� ICA_-2 ���Ϊ����(L/D-1)��Ŀ���
if ICA_down(4,1)==ECA_down(4,1) && ICA_down(4,1)~=0    %CCA1
    index(2:4)=ICA_down(4:6,2);
    r(2:4)=ICA_down(4:6,1);
elseif ICA_down(4,1)==0||ECA_down(4,1)==0
    r(2:4)=0;    
else
    cc=0.5*(data(ICA_down(4,2),1:3)+data(ECA_down(4,2),1:3));
    N=indexs(2)-index(1);
    C=repmat(cc,N,1);
%     plot3(cc(1),cc(2),cc(3),'r*');
    [dismin,index1]=min(sqrt(sum((C-data((index(1)+1:indexs(2)),1:3)).^2,2)));
    index(2)=index1+index(1);
    r(2)=data(index(2),4);
    if ICA_down(5,1)==ECA_down(5,1) && ICA_down(5,1)~=0   %CCA2
        index(3:4)=ICA_down(5:6,2);
        r(3:4)=ICA_down(5:6,1);
    elseif ICA_down(5,1)==0||ECA_down(5,1)==0
        r(3:4)=0;
    else
        cc=0.5*(data(ICA_down(5,2),1:3)+data(ECA_down(5,2),1:3));
        N=indexs(2)-index(2);
        C=repmat(cc,N,1);
%         plot3(cc(1),cc(2),cc(3),'r*');
        [dismin,index1]=min(sqrt(sum((C-data((index(2)+1:indexs(2)),1:3)).^2,2)));
        index(3)=index1+index(2);
        r(3)=data(index(3),4);
        if ICA_down(6,1)==ECA_down(6,1) &&  ICA_down(6,1)~=0  %CCA3
            index(4)=ICA_down(6,2);
            r(4)=ICA_down(6,1);
        elseif ICA_down(6,1)==0||ECA_down(6,1)==0
            r(4)=0;
        else
            cc=0.5*(data(ICA_down(6,2),1:3)+data(ECA_down(6,2),1:3));
            N=indexs(2)-index(3);
            C=repmat(cc,N,1);
%             plot3(cc(1),cc(2),cc(3),'r*');
            [dismin,index1]=min(sqrt(sum((C-data((index(3)+1:indexs(2)),1:3)).^2,2)));
            index(4)=index1+index(3);
            r(4)=data(index(4),4);
        end
    end
end
CCA=[r(1:4) index(1:4)];
%%%%%%%%%%%%%%%%% ������ݼ��� %%%%%%%%%%%
%%%%%% �Ѿ��õ� ICA/ECA/CCA/ICA_down/ECA_down
% D=sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
% ������������D = ICA1to4dis
% L=sum((sqrt(sum((data(ICA(4,2):ICA(1,2)-1,1:3)-data(ICA(4,2)+1:ICA(1,2),1:3)).^2,2))));
% ������������L = ICA14len
% ICA3_t=L/D-1;
%%%%%%%% �����жϴ��� %%%%%%%%
%%%%%% ICA��֮�䳤�Ⱥ;������ȡ&&ICA��CCA֮�䳤���������ȡ %%%%%%%%
if ICA(6,1)~=0
    ICA0to5dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
    ICA0to4dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
    ICA0to3dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
    ICA0to2dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
    ICA0to1dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
    
    ICA05len = sum((sqrt(sum((data(ICA(6,2):ICA(1,2)-1,1:3)-data(ICA(6,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA04len = sum((sqrt(sum((data(ICA(5,2):ICA(1,2)-1,1:3)-data(ICA(5,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA03len = sum((sqrt(sum((data(ICA(4,2):ICA(1,2)-1,1:3)-data(ICA(4,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA02len = sum((sqrt(sum((data(ICA(3,2):ICA(1,2)-1,1:3)-data(ICA(3,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA01len = sum((sqrt(sum((data(ICA(2,2):ICA(1,2)-1,1:3)-data(ICA(2,2)+1:ICA(1,2),1:3)).^2,2))));
    
    if CCA(4,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
        CCA2toICA5dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
        CCA1toICA5dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
        
        CCA3_ICA5len = sum((sqrt(sum((data(ICA(6,2):CCA(4,2)-1,1:3)-data(ICA(6,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA5len = sum((sqrt(sum((data(ICA(6,2):CCA(3,2)-1,1:3)-data(ICA(6,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA5len = sum((sqrt(sum((data(ICA(6,2):CCA(2,2)-1,1:3)-data(ICA(6,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        CCA2toICA4dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        CCA1toICA4dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        
        CCA3_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(4,2)-1,1:3)-data(ICA(5,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(3,2)-1,1:3)-data(ICA(5,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(2,2)-1,1:3)-data(ICA(5,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA2toICA3dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(4,2)-1,1:3)-data(ICA(4,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(3,2)-1,1:3)-data(ICA(4,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(4,2)-1,1:3)-data(ICA(3,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(4,2)-1,1:3)-data(ICA(2,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    elseif CCA(3,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;
        CCA2toICA5dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
        CCA1toICA5dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
        
        CCA3_ICA5len = 0;
        CCA2_ICA5len = sum((sqrt(sum((data(ICA(6,2):CCA(3,2)-1,1:3)-data(ICA(6,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA5len = sum((sqrt(sum((data(ICA(6,2):CCA(2,2)-1,1:3)-data(ICA(6,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;
        CCA2toICA4dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        CCA1toICA4dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        
        CCA3_ICA4len = 0;
        CCA2_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(3,2)-1,1:3)-data(ICA(5,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(2,2)-1,1:3)-data(ICA(5,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;
        CCA2toICA3dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = 0;
        CCA2_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(3,2)-1,1:3)-data(ICA(4,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
        
    elseif CCA(2,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;
        CCA1toICA5dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(6,2),1:3)).^2,2))));
        
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;
        CCA1_ICA5len = sum((sqrt(sum((data(ICA(6,2):CCA(2,2)-1,1:3)-data(ICA(6,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;
        CCA1toICA4dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;
        CCA1_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(2,2)-1,1:3)-data(ICA(5,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    else
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0; 		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;		CCA1toICA2dis = 0;
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;		CCA1_ICA2len = 0;
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;		CCA1toICA1dis = 0;
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;		CCA1_ICA1len = 0;
    end
    
elseif ICA(5,1)~=0
    ICA0to5dis = 0;
    ICA0to4dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
    ICA0to3dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
    ICA0to2dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
    ICA0to1dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
    
    ICA05len = 0;
    ICA04len = sum((sqrt(sum((data(ICA(5,2):ICA(1,2)-1,1:3)-data(ICA(5,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA03len = sum((sqrt(sum((data(ICA(4,2):ICA(1,2)-1,1:3)-data(ICA(4,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA02len = sum((sqrt(sum((data(ICA(3,2):ICA(1,2)-1,1:3)-data(ICA(3,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA01len = sum((sqrt(sum((data(ICA(2,2):ICA(1,2)-1,1:3)-data(ICA(2,2)+1:ICA(1,2),1:3)).^2,2))));
    
    if CCA(4,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        CCA2toICA4dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        CCA1toICA4dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        
        CCA3_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(4,2)-1,1:3)-data(ICA(5,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(3,2)-1,1:3)-data(ICA(5,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(2,2)-1,1:3)-data(ICA(5,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA2toICA3dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(4,2)-1,1:3)-data(ICA(4,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(3,2)-1,1:3)-data(ICA(4,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(4,2)-1,1:3)-data(ICA(3,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(4,2)-1,1:3)-data(ICA(2,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    elseif CCA(3,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;
        CCA2toICA4dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        CCA1toICA4dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        
        CCA3_ICA4len = 0;
        CCA2_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(3,2)-1,1:3)-data(ICA(5,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(2,2)-1,1:3)-data(ICA(5,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;
        CCA2toICA3dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = 0;
        CCA2_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(3,2)-1,1:3)-data(ICA(4,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
        
    elseif CCA(2,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;
        CCA1toICA4dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(5,2),1:3)).^2,2))));
        
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;
        CCA1_ICA4len = sum((sqrt(sum((data(ICA(5,2):CCA(2,2)-1,1:3)-data(ICA(5,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    else
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0; 		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;		CCA1toICA2dis = 0;
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;		CCA1_ICA2len = 0;
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;		CCA1toICA1dis = 0;
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;		CCA1_ICA1len = 0;
    end
    
elseif ICA(4,1)~=0
    ICA0to5dis = 0;
    ICA0to4dis = 0;
    ICA0to3dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
    ICA0to2dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
    ICA0to1dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
    
    ICA05len = 0;
    ICA04len = 0;
    ICA03len = sum((sqrt(sum((data(ICA(4,2):ICA(1,2)-1,1:3)-data(ICA(4,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA02len = sum((sqrt(sum((data(ICA(3,2):ICA(1,2)-1,1:3)-data(ICA(3,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA01len = sum((sqrt(sum((data(ICA(2,2):ICA(1,2)-1,1:3)-data(ICA(2,2)+1:ICA(1,2),1:3)).^2,2))));
    
    if CCA(4,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA2toICA3dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(4,2)-1,1:3)-data(ICA(4,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(3,2)-1,1:3)-data(ICA(4,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(4,2)-1,1:3)-data(ICA(3,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(4,2)-1,1:3)-data(ICA(2,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    elseif CCA(3,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;
        CCA2toICA3dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = 0;
        CCA2_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(3,2)-1,1:3)-data(ICA(4,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
        
    elseif CCA(2,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;
        CCA1toICA3dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(4,2),1:3)).^2,2))));
        
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;
        CCA1_ICA3len = sum((sqrt(sum((data(ICA(4,2):CCA(2,2)-1,1:3)-data(ICA(4,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    else
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0; 		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;		CCA1toICA2dis = 0;
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;		CCA1_ICA2len = 0;
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;		CCA1toICA1dis = 0;
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;		CCA1_ICA1len = 0;
    end
    
elseif ICA(3,1)~=0
    ICA0to5dis = 0;
    ICA0to4dis = 0;
    ICA0to3dis = 0;
    ICA0to2dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
    ICA0to1dis = sum((sqrt(sum((data(ICA(1,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
    
    ICA05len = 0;
    ICA04len = 0;
    ICA03len = 0;
    ICA02len = sum((sqrt(sum((data(ICA(3,2):ICA(1,2)-1,1:3)-data(ICA(3,2)+1:ICA(1,2),1:3)).^2,2))));
    ICA01len = sum((sqrt(sum((data(ICA(2,2):ICA(1,2)-1,1:3)-data(ICA(2,2)+1:ICA(1,2),1:3)).^2,2))));
    
    if CCA(4,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(4,2)-1,1:3)-data(ICA(3,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = sum((sqrt(sum((data(CCA(4,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(4,2)-1,1:3)-data(ICA(2,2)+1:CCA(4,2),1:3)).^2,2))));
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    elseif CCA(3,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;
        CCA2toICA2dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;
        CCA2_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(3,2)-1,1:3)-data(ICA(3,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;
        CCA2toICA1dis = sum((sqrt(sum((data(CCA(3,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;
        CCA2_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(3,2)-1,1:3)-data(ICA(2,2)+1:CCA(3,2),1:3)).^2,2))));
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
        
    elseif CCA(2,1)~=0
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0;		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;
        CCA1toICA2dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(3,2),1:3)).^2,2))));
        
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;
        CCA1_ICA2len = sum((sqrt(sum((data(ICA(3,2):CCA(2,2)-1,1:3)-data(ICA(3,2)+1:CCA(2,2),1:3)).^2,2))));
        
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;
        CCA1toICA1dis = sum((sqrt(sum((data(CCA(2,2),1:3)-data(ICA(2,2),1:3)).^2,2))));
        
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;
        CCA1_ICA1len = sum((sqrt(sum((data(ICA(2,2):CCA(2,2)-1,1:3)-data(ICA(2,2)+1:CCA(2,2),1:3)).^2,2))));
    else
        %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA5dis = 0;		CCA2toICA5dis = 0; 		CCA1toICA5dis = 0;
        CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
        %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
        CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
        %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
        CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
        %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA2dis = 0;		CCA2toICA2dis = 0;		CCA1toICA2dis = 0;
        CCA3_ICA2len = 0;		CCA2_ICA2len = 0;		CCA1_ICA2len = 0;
        %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
        CCA3toICA1dis = 0;		CCA2toICA1dis = 0;		CCA1toICA1dis = 0;
        CCA3_ICA1len = 0;		CCA2_ICA1len = 0;		CCA1_ICA1len = 0;
    end
    
else
    ICA0to5dis = 0; ICA0to4dis = 0; ICA0to3dis = 0; ICA0to2dis = 0; ICA0to1dis = 0;
    ICA05len = 0; ICA04len = 0; ICA03len = 0; ICA02len = 0; ICA01len = 0;
    %%%%%% CCAxx_ICA5_���Ⱥ;������ȡ %%%%%%%%
    CCA3toICA5dis = 0;		CCA2toICA5dis = 0; 		CCA1toICA5dis = 0;
    CCA3_ICA5len = 0;		CCA2_ICA5len = 0;		CCA1_ICA5len = 0;
    %%%%%% CCAxx_ICA4_���Ⱥ;������ȡ %%%%%%%%
    CCA3toICA4dis = 0;		CCA2toICA4dis = 0;		CCA1toICA4dis = 0;
    CCA3_ICA4len = 0;		CCA2_ICA4len = 0;		CCA1_ICA4len = 0;
    %%%%%% CCAxx_ICA3_���Ⱥ;������ȡ %%%%%%%%
    CCA3toICA3dis = 0;		CCA2toICA3dis = 0;		CCA1toICA3dis = 0;
    CCA3_ICA3len = 0;		CCA2_ICA3len = 0;		CCA1_ICA3len = 0;
    %%%%%% CCAxx_ICA2_���Ⱥ;������ȡ %%%%%%%%
    CCA3toICA2dis = 0;		CCA2toICA2dis = 0;		CCA1toICA2dis = 0;
    CCA3_ICA2len = 0;		CCA2_ICA2len = 0;		CCA1_ICA2len = 0;
    %%%%%% CCAxx_ICA1_���Ⱥ;������ȡ %%%%%%%%
    CCA3toICA1dis = 0;		CCA2toICA1dis = 0;		CCA1toICA1dis = 0;
    CCA3_ICA1len = 0;		CCA2_ICA1len = 0;		CCA1_ICA1len = 0;
end

%%%%%%%%%%%%%%
if CCA(4,1)~=0
    %%%%%% CCA��֮�䳤�Ⱥ;������ȡ %%%%%%%%
    CCA0to3dis = sum((sqrt(sum((data(CCA(1,2),1:3)-data(CCA(4,2),1:3)).^2,2))));
    CCA0to2dis = sum((sqrt(sum((data(CCA(1,2),1:3)-data(CCA(3,2),1:3)).^2,2))));
    CCA0to1dis = sum((sqrt(sum((data(CCA(1,2),1:3)-data(CCA(2,2),1:3)).^2,2))));
    
    CCA03len = sum((sqrt(sum((data(CCA(1,2):CCA(4,2)-1,1:3)-data(CCA(1,2)+1:CCA(4,2),1:3)).^2,2))));
    CCA02len = sum((sqrt(sum((data(CCA(1,2):CCA(3,2)-1,1:3)-data(CCA(1,2)+1:CCA(3,2),1:3)).^2,2))));
    CCA01len = sum((sqrt(sum((data(CCA(1,2):CCA(2,2)-1,1:3)-data(CCA(1,2)+1:CCA(2,2),1:3)).^2,2))));
    
    %%%%%% CCAxx_ECA1_���Ⱥ;������ȡ %%%%%%%% ���ο�!!!!
    CCA3toECA1dis = sum((sqrt(sum((data(ECA_down(6,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
    CCA2toECA1dis = sum((sqrt(sum((data(ECA_down(5,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
    CCA1toECA1dis = sum((sqrt(sum((data(ECA_down(4,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
    
    CCA3_ECA1len = sum((sqrt(sum((data(ECA(2,2):ECA_down(6,2)-1,1:3)-data(ECA(2,2)+1:ECA_down(6,2),1:3)).^2,2))));
    CCA2_ECA1len = sum((sqrt(sum((data(ECA(2,2):ECA_down(5,2)-1,1:3)-data(ECA(2,2)+1:ECA_down(5,2),1:3)).^2,2))));
    CCA1_ECA1len = sum((sqrt(sum((data(ECA(2,2):ECA_down(4,2)-1,1:3)-data(ECA(2,2)+1:ECA_down(4,2),1:3)).^2,2))));
elseif CCA(3,1)~=0
    %%%%%% CCA��֮�䳤�Ⱥ;������ȡ %%%%%%%%
    CCA0to3dis = 0;
    CCA0to2dis = sum((sqrt(sum((data(CCA(1,2),1:3)-data(CCA(3,2),1:3)).^2,2))));
    CCA0to1dis = sum((sqrt(sum((data(CCA(1,2),1:3)-data(CCA(2,2),1:3)).^2,2))));
    
    CCA03len = 0;
    CCA02len = sum((sqrt(sum((data(CCA(1,2):CCA(3,2)-1,1:3)-data(CCA(1,2)+1:CCA(3,2),1:3)).^2,2))));
    CCA01len = sum((sqrt(sum((data(CCA(1,2):CCA(2,2)-1,1:3)-data(CCA(1,2)+1:CCA(2,2),1:3)).^2,2))));
    
    %%%%%% CCAxx_ECA1_���Ⱥ;������ȡ %%%%%%%% ���ο�!!!!
    CCA3toECA1dis = 0;
    CCA2toECA1dis = sum((sqrt(sum((data(ECA_down(5,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
    CCA1toECA1dis = sum((sqrt(sum((data(ECA_down(4,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
    CCA3_ECA1len = 0;
    CCA2_ECA1len = sum((sqrt(sum((data(ECA(2,2):ECA_down(5,2)-1,1:3)-data(ECA(2,2)+1:ECA_down(5,2),1:3)).^2,2))));
    CCA1_ECA1len = sum((sqrt(sum((data(ECA(2,2):ECA_down(4,2)-1,1:3)-data(ECA(2,2)+1:ECA_down(4,2),1:3)).^2,2))));
elseif CCA(2,1)~=0
    %%%%%% CCA��֮�䳤�Ⱥ;������ȡ %%%%%%%%
    CCA0to3dis = 0;		CCA0to2dis = 0;
    CCA0to1dis = sum((sqrt(sum((data(CCA(1,2),1:3)-data(CCA(2,2),1:3)).^2,2))));
    
    CCA03len = 0;		CCA02len = 0;
    CCA01len = sum((sqrt(sum((data(CCA(1,2):CCA(2,2)-1,1:3)-data(CCA(1,2)+1:CCA(2,2),1:3)).^2,2))));
    
    %%%%%% CCAxx_ECA1_���Ⱥ;������ȡ %%%%%%%% ���ο�!!!!
    CCA3toECA1dis = 0;		CCA2toECA1dis = 0;
    CCA1toECA1dis = sum((sqrt(sum((data(ECA_down(4,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
    CCA3_ECA1len = 0;		CCA2_ECA1len = 0;
    CCA1_ECA1len = sum((sqrt(sum((data(ECA(2,2):ECA_down(4,2)-1,1:3)-data(ECA(2,2)+1:ECA_down(4,2),1:3)).^2,2))));
else
    %%%%%% CCA��֮�䳤�Ⱥ;������ȡ %%%%%%%%
    CCA0to3dis = 0;		CCA0to2dis = 0;		CCA0to1dis = 0;
    CCA03len = 0; 		CCA02len = 0;		CCA01len = 0;
    %%%%%% CCAxx_ECA1_���Ⱥ;������ȡ %%%%%%%% ���ο�!!!!
    CCA3toECA1dis = 0;		CCA2toECA1dis = 0;		CCA1toECA1dis = 0;
    CCA3_ECA1len = 0;		CCA2_ECA1len = 0;		CCA1_ECA1len = 0;
end
%%%%%% ECA��֮�䳤�Ⱥ;������ȡ %%%%%%%%
ECA0to1dis = sum((sqrt(sum((data(ECA(1,2),1:3)-data(ECA(2,2),1:3)).^2,2))));
ECA01len = sum((sqrt(sum((data(ECA(2,2):ECA(1,2)-1,1:3)-data(ECA(2,2)+1:ECA(1,2),1:3)).^2,2))));

%�������
results={'ICA0to5dis', 'ICA0to4dis','ICA0to3dis',   'ICA0to2dis',   'ICA0to1dis',...
         'ICA05len',   'ICA04len',  'ICA03len',     'ICA02len',     'ICA01len',...
         'ECA0to1dis', 'ECA01len',  'CCA3toECA1dis','CCA2toECA1dis','CCA1toECA1dis',...
         'CCA3_ECA1len','CCA2_ECA1len','CCA1_ECA1len',...
         'CCA0to3dis','CCA0to2dis','CCA0to1dis', 'CCA03len','CCA02len','CCA01len',...
         'CCA3toICA5dis','CCA2toICA5dis','CCA1toICA5dis','CCA3_ICA5len','CCA2_ICA5len','CCA1_ICA5len',...
         'CCA3toICA4dis','CCA2toICA4dis','CCA1toICA4dis','CCA3_ICA4len','CCA2_ICA4len','CCA1_ICA4len',...
         'CCA3toICA3dis','CCA2toICA3dis','CCA1toICA3dis','CCA3_ICA3len','CCA2_ICA3len','CCA1_ICA3len',...
         'CCA3toICA2dis','CCA2toICA2dis','CCA1toICA2dis','CCA3_ICA2len','CCA2_ICA2len','CCA1_ICA2len',...
         'CCA3toICA1dis','CCA2toICA1dis','CCA1toICA1dis','CCA3_ICA1len','CCA2_ICA1len','CCA1_ICA1len';...
         ICA0to5dis,ICA0to4dis,ICA0to3dis,ICA0to2dis,ICA0to1dis,...
         ICA05len,  ICA04len,  ICA03len,  ICA02len,  ICA01len ,...
         ECA0to1dis,ECA01len,  CCA3toECA1dis,CCA2toECA1dis,CCA1toECA1dis ,...
         CCA3_ECA1len,CCA2_ECA1len,CCA1_ECA1len,...
         CCA0to3dis, CCA0to2dis, CCA0to1dis, CCA03len, CCA02len, CCA01len,...
         CCA3toICA5dis,CCA2toICA5dis,CCA1toICA5dis,CCA3_ICA5len,CCA2_ICA5len,CCA1_ICA5len,...
         CCA3toICA4dis,CCA2toICA4dis,CCA1toICA4dis,CCA3_ICA4len,CCA2_ICA4len,CCA1_ICA4len,...
         CCA3toICA3dis,CCA2toICA3dis,CCA1toICA3dis,CCA3_ICA3len,CCA2_ICA3len,CCA1_ICA3len,...
         CCA3toICA2dis,CCA2toICA2dis,CCA1toICA2dis,CCA3_ICA2len,CCA2_ICA2len,CCA1_ICA2len,...
         CCA3toICA1dis,CCA2toICA1dis,CCA1toICA1dis,CCA3_ICA1len,CCA2_ICA1len,CCA1_ICA1len }; % 2-8cell
pointradius ={ 'ICA5radius','ICA4radius','ICA3radius','ICA2radius','ICA1radius',...
    'ECA1radius','CCA3radius','CCA2radius','CCA1radius';
    ICA(6,1),ICA(5,1),ICA(4,1),ICA(3,1),ICA(2,1),ECA(2,1),CCA(4,1),CCA(3,1),CCA(2,1)};
end