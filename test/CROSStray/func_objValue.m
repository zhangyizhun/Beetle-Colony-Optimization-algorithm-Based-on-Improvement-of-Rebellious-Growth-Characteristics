function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明

fact1 = sin(x).*sin(y);
fact2 = exp(abs(100 - sqrt(x.^2+y.^2)./pi));

z = -0.0001 .* (abs(fact1.*fact2)+1).^0.1;

result = z ;
 
end