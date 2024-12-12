function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明
d = 2;


w1 = 1 + (x - 1)./4;
w2=1 + (y - 1)./4;

term1 = (sin(pi.*w1)).^2;
term3 = (w2-1).^2 .* (1+(sin(2.*pi.*w2)).^2);

sum = 0;

        new = (w1-1).^2 .* (1+10.*(sin(pi.*w1+1)).^2);
	sum = sum + new;


z = term1 + sum + term3;

result = z ;
 
end