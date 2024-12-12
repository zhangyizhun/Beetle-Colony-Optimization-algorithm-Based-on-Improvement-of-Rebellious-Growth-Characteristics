function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明

temp1=x.^2+y.^2;

temp2=cos(2*pi*x)+cos(2*pi*y);

z=20+exp(1)-20*exp(-0.2*sqrt(temp1/2))-exp(temp2/2);

result = z ;
 
end