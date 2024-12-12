function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明

z=x.^2+y.^2-10*cos(2*pi*x)-10*cos(2*pi*y)+20;

result = z ;
 
end