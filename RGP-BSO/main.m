clc 
clear all
%% I.绘制目标函数曲线

[x,y]=meshgrid(-5:0.1:5,-5:0.1:5);
z=x.^2+y.^2-10*cos(2*pi*x)-10*cos(2*pi*y)+20;
mesh(x,y,z);
hold on
%% II.参数初始化

c1=1.49445;
c2=1.49445;

 ws=0.9;
 we=0.4;
 k=0.4;%比例

maxgen=20;%进化次数
sizepop=50;%种群规模
Vmax=1;%速度最大值和速度最小值
Vmin=-1;

popmax=5;%飞行范围根据函数定义到1-2
popmin=-5;

 %初始步长
 step=1;
c=1;
 %% 产生初始粒子和速度
for i=1:sizepop
%随机产生种群
pop(i,:)=5*rands(1,2);%定义到-5,5  rans中是后边为2因为是二元函数
V(i,:)=rands(1,2);%%初始化速度，与定义的最大值最小值要一致
%%%%计算适应度
fitness(i)=fun(pop(i,1),pop(i,2));
end

%% 个体极值和群体极值
[bestfitness bestindex]=max(fitness);
gbest=pop(bestindex,:);%全局最佳
pbest=pop;%个体最佳
fitnesspbest = fitness;%个体最佳适应度值；
fitnessgbest = bestfitness; %全局最佳适应度值


%% 迭代寻优

for i=1:maxgen
    w=ws-(ws-we)*(i/maxgen);
  for j=1:sizepop
    %天牛子群部分位置移动
       d0=step/c;%两须之间的距离
        xleft=pop(j,:)+V(j,:)*d0/2;
        fleft=fun(xleft(1,1),xleft(1,2));
        xright=pop(j,:)-V(j,:)*d0/2;
        fright=fun(xright(1,1),xright(1,2));
        Y(j,:)=step.*V(j,:).*sign(fleft-fright);


   % V(j,:)=V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%速度更新公式
     V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%速度更新公式
      


      V(j,find(V(j,:)>Vmax))=Vmax;%速度约束
      V(j,find(V(j,:)<Vmin))=Vmin;

      % 种群更新
      pop(j,:)=pop(j,:)+k*V(j,:)+(1-k)*Y(j,:);

%叛逆性格
 k=k+0.006;
    if k>1
        k=1;
    end

      pop(j,find(pop(j,:)>popmax))=popmax;%位置约束
      pop(j,find(pop(j,:)<popmin))=popmin;

      %适应度值更新
      fitness(j)=fun(pop(j,1),pop(j,2));
  end

  for j=1:sizepop
     %个体最优值更新
     if fitness(j)>fitnesspbest(j)
          pbest(j,:)=pop(j,:);
          fitnesspbest(j)=fitness(j);
     end
  %群体最优值更新
  if fitness(j)>fitnessgbest
      gbest=pop(j,:);
      fitnessgbest=fitness(j);
     
  end
     
  end
yy(i)=fitnessgbest;
end


%% 输出结果并绘图
[fitnessgbest gbest]
plot3(gbest(1),gbest(2),fitnessgbest,'bo','linewidth',1.5);
figure
plot(yy)
title('最优个体适应值','fontsize',12);
xlabel('进化代数','fontsize',12);
ylabel('适应度','fontsize',12);

