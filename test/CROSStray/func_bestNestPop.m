function [ nestPop ] = func_bestNestPop( nestPop,new_nestPop )
%FUNC_ 此处显示有关此函数的摘要
%@author zhaoyuqiang
index = find(func_fitness(nestPop)>func_fitness(new_nestPop)) ;
nestPop(index,:) = new_nestPop(index,:) ;
end
