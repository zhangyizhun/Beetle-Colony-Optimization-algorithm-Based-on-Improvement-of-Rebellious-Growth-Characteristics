function z = fun(x,y)
%FUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


z=x.^2+y.^2-10*cos(2*pi*x)-10*cos(2*pi*y)+20;

end
% 
% function yout=fun(x)
% 
% theta=x;
% x=theta(1);
% y=theta(2);
% %yout=-sin(x).*(sin(x.^2/pi)).^20-sin(y).*(sin(2*y.^2/pi)).^20;
% %yout=x.^2+y.^2-10*cos(2*pi*x)-10*cos(2*pi*y)+20;
% 
% yout=sin(10*pi*x)/x
% end
