function z = fun(x,y)
%FUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


fact1 = sin(x).*cos(y);
fact2 = exp(abs(1 - sqrt(x.^2+y.^2)./pi));

z = -abs(fact1.*fact2);
end

