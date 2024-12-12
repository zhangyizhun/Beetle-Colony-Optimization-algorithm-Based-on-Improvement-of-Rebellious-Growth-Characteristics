function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明


frac1 = 1 + cos(12.*sqrt(x.^2+y.^2));
frac2 = 0.5.*(x.^2+y.^2) + 2;

z = -frac1./frac2;

result = z ;
 
end