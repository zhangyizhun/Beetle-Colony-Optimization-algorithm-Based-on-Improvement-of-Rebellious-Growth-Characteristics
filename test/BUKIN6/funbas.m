function z = funbas(x,y)
%FUN 此处显示有关此函数的摘要
%   此处显示详细说明
term1 = 100 * sqrt(abs(y - 0.01*x^2));
term2 = 0.01 * abs(x+10);

z = term1 + term2;

z=-z;

end

