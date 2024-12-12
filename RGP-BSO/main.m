clc 
clear all
%% I.����Ŀ�꺯������

[x,y]=meshgrid(-5:0.1:5,-5:0.1:5);
z=x.^2+y.^2-10*cos(2*pi*x)-10*cos(2*pi*y)+20;
mesh(x,y,z);
hold on
%% II.������ʼ��

c1=1.49445;
c2=1.49445;

 ws=0.9;
 we=0.4;
 k=0.4;%����

maxgen=20;%��������
sizepop=50;%��Ⱥ��ģ
Vmax=1;%�ٶ����ֵ���ٶ���Сֵ
Vmin=-1;

popmax=5;%���з�Χ���ݺ������嵽1-2
popmin=-5;

 %��ʼ����
 step=1;
c=1;
 %% ������ʼ���Ӻ��ٶ�
for i=1:sizepop
%���������Ⱥ
pop(i,:)=5*rands(1,2);%���嵽-5,5  rans���Ǻ��Ϊ2��Ϊ�Ƕ�Ԫ����
V(i,:)=rands(1,2);%%��ʼ���ٶȣ��붨������ֵ��СֵҪһ��
%%%%������Ӧ��
fitness(i)=fun(pop(i,1),pop(i,2));
end

%% ���弫ֵ��Ⱥ�弫ֵ
[bestfitness bestindex]=max(fitness);
gbest=pop(bestindex,:);%ȫ�����
pbest=pop;%�������
fitnesspbest = fitness;%���������Ӧ��ֵ��
fitnessgbest = bestfitness; %ȫ�������Ӧ��ֵ


%% ����Ѱ��

for i=1:maxgen
    w=ws-(ws-we)*(i/maxgen);
  for j=1:sizepop
    %��ţ��Ⱥ����λ���ƶ�
       d0=step/c;%����֮��ľ���
        xleft=pop(j,:)+V(j,:)*d0/2;
        fleft=fun(xleft(1,1),xleft(1,2));
        xright=pop(j,:)-V(j,:)*d0/2;
        fright=fun(xright(1,1),xright(1,2));
        Y(j,:)=step.*V(j,:).*sign(fleft-fright);


   % V(j,:)=V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%�ٶȸ��¹�ʽ
     V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%�ٶȸ��¹�ʽ
      


      V(j,find(V(j,:)>Vmax))=Vmax;%�ٶ�Լ��
      V(j,find(V(j,:)<Vmin))=Vmin;

      % ��Ⱥ����
      pop(j,:)=pop(j,:)+k*V(j,:)+(1-k)*Y(j,:);

%�����Ը�
 k=k+0.006;
    if k>1
        k=1;
    end

      pop(j,find(pop(j,:)>popmax))=popmax;%λ��Լ��
      pop(j,find(pop(j,:)<popmin))=popmin;

      %��Ӧ��ֵ����
      fitness(j)=fun(pop(j,1),pop(j,2));
  end

  for j=1:sizepop
     %��������ֵ����
     if fitness(j)>fitnesspbest(j)
          pbest(j,:)=pop(j,:);
          fitnesspbest(j)=fitness(j);
     end
  %Ⱥ������ֵ����
  if fitness(j)>fitnessgbest
      gbest=pop(j,:);
      fitnessgbest=fitness(j);
     
  end
     
  end
yy(i)=fitnessgbest;
end


%% ����������ͼ
[fitnessgbest gbest]
plot3(gbest(1),gbest(2),fitnessgbest,'bo','linewidth',1.5);
figure
plot(yy)
title('���Ÿ�����Ӧֵ','fontsize',12);
xlabel('��������','fontsize',12);
ylabel('��Ӧ��','fontsize',12);

