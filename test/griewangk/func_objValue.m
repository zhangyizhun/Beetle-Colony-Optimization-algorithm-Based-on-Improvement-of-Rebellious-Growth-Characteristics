function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明

[x,y]=meshgrid(x);

temp1=(x.^2+y.^2)/4000;

temp2=cos(x)*cos(y/sqrt(2));

z=temp1-temp2+1;

result = z ;
 
end