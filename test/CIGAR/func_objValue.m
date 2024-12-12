function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明

z=x.^2+(10^4)*y.^2;

result = z ;
 
end