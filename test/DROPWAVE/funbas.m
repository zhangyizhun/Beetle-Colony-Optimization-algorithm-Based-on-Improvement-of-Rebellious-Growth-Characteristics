function z = funbas(x,y)
%FUN 此处显示有关此函数的摘要
%   此处显示详细说明



frac1 = 1 + cos(12.*sqrt(x.^2+y.^2));
frac2 = 0.5.*(x.^2+y.^2) + 2;

z = -frac1./frac2;
z=-z;

end

