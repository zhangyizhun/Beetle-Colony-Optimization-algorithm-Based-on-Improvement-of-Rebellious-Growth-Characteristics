function z = fun(x,y)
%FUN 此处显示有关此函数的摘要
%   此处显示详细说明


temp1=x.^2+y.^2;

temp2=cos(2*pi*x)+cos(2*pi*y);

z=20+exp(1)-20*exp(-0.2*sqrt(temp1/2))-exp(temp2/2);

end

