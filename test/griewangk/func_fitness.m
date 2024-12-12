function [ result ] = func_fitness(x)
%OBJFUNCTION 求适应度，最小值
%   待优化目标函数
% x:　种群或者个体
% result : 种群适应度
 theta=x;
 x=theta(1);
 y=theta(2);
objValue =  func_objValue(x,y);
result  = objValue ;
end