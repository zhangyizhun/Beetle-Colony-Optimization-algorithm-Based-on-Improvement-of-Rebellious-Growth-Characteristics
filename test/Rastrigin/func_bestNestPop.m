function [ nestPop ] = func_bestNestPop( nestPop,new_nestPop )
%FUNC_ �˴���ʾ�йش˺�����ժҪ
%@author zhaoyuqiang
index = find(func_fitness(nestPop)>func_fitness(new_nestPop)) ;
nestPop(index,:) = new_nestPop(index,:) ;
end
