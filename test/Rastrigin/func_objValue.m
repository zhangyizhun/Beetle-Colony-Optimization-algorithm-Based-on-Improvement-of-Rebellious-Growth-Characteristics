function [ result ] = func_objValue(x,y )
%FUNC_OBJVALUE ����Ŀ�꺯��
%   �˴���ʾ��ϸ˵��

z=x.^2+y.^2-10*cos(2*pi*x)-10*cos(2*pi*y)+20;

result = z ;
 
end