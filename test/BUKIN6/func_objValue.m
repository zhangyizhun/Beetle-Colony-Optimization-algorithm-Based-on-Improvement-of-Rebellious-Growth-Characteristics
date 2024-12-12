function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明
term1 = 100 * sqrt(abs(y - 0.01*x^2));
term2 = 0.01 * abs(x+10);

z = term1 + term2;


result = z ;
 
end