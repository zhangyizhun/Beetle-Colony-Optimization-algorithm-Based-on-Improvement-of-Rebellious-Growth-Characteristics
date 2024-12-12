function [ result ] = func_levy( nestPop,Xmax,Xmin)
%FUNC_LEVY : Update position of nest by using Levy flights
%@author : zhaoyuqiang 
[N,D] = size(nestPop) ;
% Levy flights by Mantegna's algorithm
beta = 1.5 ;
alpha = 1 ;
sigma_u = (gamma(1+beta)*sin(pi*beta/2)/(beta*gamma((1+beta)/2)*2^((beta-1)/2)))^(1/beta) ;
sigma_v = 1 ;
u = normrnd(0,sigma_u,N,D) ;
v = normrnd(0,sigma_v,N,D) ;
step = u./(abs(v).^(1/beta)) ;
% alpha = 0.1.*(nestPop(randperm(N),:)-nestPop(randperm(N),:)); % Bad effect
nestPop = nestPop+alpha.*step ;
% Deal with bounds
nestPop(find(nestPop>Xmax)) = Xmax ; %#ok<*FNDSB>
nestPop(find(nestPop<Xmin)) = Xmin ;
result = nestPop ; 
end
