function z = funbas(x,y)
%FUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

[x,y]=meshgrid(x);

temp1=(x.^2+y.^2)/4000;

temp2=cos(x)*cos(y/sqrt(2));

z=temp1-temp2+1;
z=-z;

end

