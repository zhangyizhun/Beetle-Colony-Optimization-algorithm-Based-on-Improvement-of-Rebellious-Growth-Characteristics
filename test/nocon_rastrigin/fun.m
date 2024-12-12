function z = fun(x,y)
%FUN 此处显示有关此函数的摘要
%   此处显示详细说明
temp1=0;

temp2=0;

z=0;

if abs(x)<1/2

temp1=x;

else

temp1=round(2.*x)/2;

end

if abs(y)<1/2

temp2=y;

else

temp2=round(2.*y)/2;

end

z=z+temp1.^2-10*cos(2*pi.*temp1)+temp2.^2-10*cos(2*pi.*temp2)+20;

end

