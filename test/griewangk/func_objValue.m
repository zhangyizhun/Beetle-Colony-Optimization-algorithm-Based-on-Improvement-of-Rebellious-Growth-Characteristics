function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE ����Ŀ�꺯��
%   �˴���ʾ��ϸ˵��

[x,y]=meshgrid(x);

temp1=(x.^2+y.^2)/4000;

temp2=cos(x)*cos(y/sqrt(2));

z=temp1-temp2+1;

result = z ;
 
end