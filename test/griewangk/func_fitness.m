function [ result ] = func_fitness(x)
%OBJFUNCTION ����Ӧ�ȣ���Сֵ
%   ���Ż�Ŀ�꺯��
% x:����Ⱥ���߸���
% result : ��Ⱥ��Ӧ��
 theta=x;
 x=theta(1);
 y=theta(2);
objValue =  func_objValue(x,y);
result  = objValue ;
end