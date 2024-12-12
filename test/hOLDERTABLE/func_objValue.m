function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明


fact1 = sin(x).*cos(y);
fact2 = exp(abs(1 - sqrt(x.^2+y.^2)./pi));

z = -abs(fact1.*fact2);

result = z ;
 
end