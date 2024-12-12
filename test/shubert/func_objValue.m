function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明
sum1 = 0;
sum2 = 0;

for ii = 1:5
	new1 = ii .* cos((ii+1).*x+ii);
	new2 = ii .* cos((ii+1).*y+ii);
	sum1 = sum1 + new1;
	sum2 = sum2 + new2;
end

z = sum1 .* sum2;

result = z ;
 
end