function z = fun(x,y)
%FUN 此处显示有关此函数的摘要
%   此处显示详细说明


fact1 = sin(x).*cos(y);
fact2 = exp(abs(1 - sqrt(x.^2+y.^2)./pi));

z = -abs(fact1.*fact2);
end

