function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE 计算目标函数
%   此处显示详细说明
z=0;

for k=0:20

z=z+0.5^k*cos(2*pi*3^k*(x+0.5))+0.5^k*cos(2*pi*3^k*(y+0.5));

end

temp=0;

for k=0:20

temp=temp+0.5^k*cos(2*pi*3^k*0.5);

end

z=z-2*temp;

result = z ;
 
end