function z = fun(x,y)
%FUN 此处显示有关此函数的摘要
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

end

