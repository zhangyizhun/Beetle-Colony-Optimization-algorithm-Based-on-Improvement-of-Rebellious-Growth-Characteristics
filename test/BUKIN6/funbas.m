function z = funbas(x,y)
%FUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
term1 = 100 * sqrt(abs(y - 0.01*x^2));
term2 = 0.01 * abs(x+10);

z = term1 + term2;

z=-z;

end

