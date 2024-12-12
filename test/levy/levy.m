
clear
clc

%levf  ��СֵΪ0  ��1,1ʱ

x=-10:0.01:10;
y=-10:0.01:10;
[x,y]=meshgrid(x);

d = 2;


w1 = 1 + (x - 1)./4;
w2=1 + (y - 1)./4;

term1 = (sin(pi.*w1)).^2;
term3 = (w2-1).^2 .* (1+(sin(2.*pi.*w2)).^2);

sum = 0;

        new = (w1-1).^2 .* (1+10.*(sin(pi.*w1+1)).^2);
	sum = sum + new;


z = term1 + sum + term3;


axis([-10,10,-10,10]);
% 
% meshc(x,y,z);
% hold on

%% 
%% pso
%% II.������ʼ��

c1=1.49445;
c2=1.49445;

 ws=0.9;
 we=0.4;



maxgen=100;%��������
sizepop=50;%��Ⱥ��ģ
Vmax=1;%�ٶ����ֵ���ٶ���Сֵ
Vmin=-1;

popmax=10;%���з�Χ���ݺ������嵽1-2
popmin=-10;


 %% ������ʼ���Ӻ��ٶ�

for i=1:sizepop
%���������Ⱥ
pop(i,:)=10*rands(1,2);%���嵽-5,5  rans���Ǻ��Ϊ2��Ϊ�Ƕ�Ԫ����
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
     % V(j,:)=V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%�ٶȸ��¹�ʽ
     V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%�ٶȸ��¹�ʽ
      V(j,find(V(j,:)>Vmax))=Vmax;%�ٶ�Լ��
      V(j,find(V(j,:)<Vmin))=Vmin;

      % ��Ⱥ����
      pop(j,:)=pop(j,:)+V(j,:);
      pop(j,find(pop(j,:)>popmax))=popmax;%λ��Լ��
      pop(j,find(pop(j,:)<popmin))=popmin;

      %��Ӧ��ֵ����
      fitness(j)=fun(pop(j,1),pop(j,2));
  end

  for j=1:sizepop
     %��������ֵ����
     if fitness(j)<fitnesspbest(j)
          pbest(j,:)=pop(j,:);
          fitnesspbest(j)=fitness(j);
     end
  %Ⱥ������ֵ����
  if fitness(j)<fitnessgbest
      gbest=pop(j,:);
      fitnessgbest=fitness(j);
     
  end
     
  end
yy(i)=fitnessgbest;
end




%% BAS

%% II.������ʼ��
%antenna distance
d0=0.001;
d1=3;
d=d1;
eta_d=0.95;
c=1;%�������ʼ����֮��Ĺ�ϵ
%random walk
l0=0.0;
l1=0.0;
l=l1;
eta_l=0.95;
%steps
step=1;%step length��ʼ����
eta_step=0.95;
n=100;%iterations ��������
k=2;%space dimension  
x0=10*rands(k,1);
%x0(1,1)=3.522993650155604;
%x0(2,1)=4.522993659400573;
x=x0;
xbest=x0;
fbest=funbas(xbest(1,1),xbest(2,1));
fbest_store=fbest;
x_store=[0;x;fbest];
display(['0:','xbest=[',num2str(xbest(1)),num2str(xbest(2)),'],fbest=',num2str(fbest)])
%% ����Ѱ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    d0=step/c;
    dir=rands(k,1);
    dir=dir/(eps+norm(dir));
    xleft=x+dir*d0;
    fleft=funbas(xleft(1,1),xleft(2,1));
    xright=x-dir*d0;
    fright=funbas(xright(1,1),xright(2,1));
    w=l*rands(k,1);
    x=x+step*dir*sign(fleft-fright);

for z=1:k
if x(z,:)>10
x(z,:)=10;
end
if x(z,:)<-10
x(z,:)=-10;
end
end
    f=funbas(x(1,1),x(2,1));
    %%%%%%%%%%%
    if f>fbest
        xbest=x;
        fbest=f;
    end

% for i=1:k
% if xbest(i,1)<-5
% xbest(i,1)=-5;
% end
% if xbest(i,1)>5
%  xbest(i,1)=5;
% end
% 
% end
    %%%%%%%%%%%
    x_store=cat(2,x_store,[i;x;f]);
    fbest_store=[fbest_store;fbest];
    display([num2str(i),':xbest=[',num2str(xbest(1)),num2str(xbest(2)),'],fbest=',num2str(fbest)])
    %%%%%%%%%%%
    d=d*eta_d+d0;
    l=l*eta_l+l0;
    step=step*eta_step;

end


%% BSO


%% II.������ʼ��

c1=1.49445;
c2=1.49445;

 ws=0.9;
 we=0.4;
 k=0.4;%����

maxgen=100;%��������
sizepop=50;%��Ⱥ��ģ
Vmax=1;%�ٶ����ֵ���ٶ���Сֵ
Vmin=-1;

popmax=10;%���з�Χ���ݺ������嵽1-2
popmin=-10;

 %��ʼ����
 step=1;
c=1;
 %% ������ʼ���Ӻ��ٶ�
for i=1:sizepop
%���������Ⱥ
pop(i,:)=10*rands(1,2);%���嵽-5,5  rans���Ǻ��Ϊ2��Ϊ�Ƕ�Ԫ����
V(i,:)=rands(1,2);%%��ʼ���ٶȣ��붨������ֵ��СֵҪһ��
%%%%������Ӧ��
fitness(i)=funbas(pop(i,1),pop(i,2));
end

%% ���弫ֵ��Ⱥ�弫ֵ
[bestfitness bestindex]=max(fitness);
gbestbso=pop(bestindex,:);%ȫ�����
pbest=pop;%�������
fitnesspbest = fitness;%���������Ӧ��ֵ��
fitnessgbestbso = bestfitness; %ȫ�������Ӧ��ֵ


%% ����Ѱ��

for i=1:maxgen
    w=ws-(ws-we)*(i/maxgen);
  for j=1:sizepop
    %��ţ��Ⱥ����λ���ƶ�
       d0=step/c;%����֮��ľ���
        xleft=pop(j,:)+V(j,:)*d0/2;
        fleft=funbas(xleft(1,1),xleft(1,2));
        xright=pop(j,:)-V(j,:)*d0/2;
        fright=funbas(xright(1,1),xright(1,2));
        Y(j,:)=step.*V(j,:).*sign(fleft-fright);


   % V(j,:)=V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%�ٶȸ��¹�ʽ
     V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbestbso-pop(j,:));%%�ٶȸ��¹�ʽ
      


      V(j,find(V(j,:)>Vmax))=Vmax;%�ٶ�Լ��
      V(j,find(V(j,:)<Vmin))=Vmin;

      % ��Ⱥ����
      pop(j,:)=pop(j,:)+k*V(j,:)+(1-k)*Y(j,:);
      pop(j,find(pop(j,:)>popmax))=popmax;%λ��Լ��
      pop(j,find(pop(j,:)<popmin))=popmin;

      %��Ӧ��ֵ����
      fitness(j)=funbas(pop(j,1),pop(j,2));
  end

  for j=1:sizepop
     %��������ֵ����
     if fitness(j)>fitnesspbest(j)
          pbest(j,:)=pop(j,:);
          fitnesspbest(j)=fitness(j);
     end
  %Ⱥ������ֵ����
  if fitness(j)>fitnessgbestbso
      gbest=pop(j,:);
      fitnessgbestbso=fitness(j);
     
  end
     
  end
bso(i)=fitnessgbestbso;
end


%% �Ľ���BSO
%% II.������ʼ��

c1=1.49445;
c2=1.49445;

 ws=0.9;
 we=0.4;
 k=0.5;%����

maxgen=100;%��������
sizepop=50;%��Ⱥ��ģ
Vmax=1;%�ٶ����ֵ���ٶ���Сֵ
Vmin=-1;

popmax=10;%���з�Χ���ݺ������嵽1-2
popmin=-10;

 %��ʼ����
 step=1;
c=1;
 %% ������ʼ���Ӻ��ٶ�
for i=1:sizepop
%���������Ⱥ
pop(i,:)=10*rands(1,2);%���嵽-5,5  rans���Ǻ��Ϊ2��Ϊ�Ƕ�Ԫ����
V(i,:)=rands(1,2);%%��ʼ���ٶȣ��붨������ֵ��СֵҪһ��
%%%%������Ӧ��
fitness(i)=funbas(pop(i,1),pop(i,2));
end

%% ���弫ֵ��Ⱥ�弫ֵ
[bestfitness bestindex]=max(fitness);
gbestgbso=pop(bestindex,:);%ȫ�����
pbest=pop;%�������
fitnesspbest = fitness;%���������Ӧ��ֵ��
fitnessgbestgbso = bestfitness; %ȫ�������Ӧ��ֵ


%% ����Ѱ��

id=randperm(sizepop);
for i=1:maxgen
    w=ws-(ws-we)*(i/maxgen);
  for j=1:sizepop
    %��ţ��Ⱥ����λ���ƶ�
       d0=step/c;%����֮��ľ���
        xleft=pop(id(j),:)+V(j,:)*d0/2;
        fleft=funbas(xleft(1,1),xleft(1,2));
        xright=pop(id(j),:)-V(j,:)*d0/2;
        fright=funbas(xright(1,1),xright(1,2));
        Y(j,:)=step.*V(j,:).*sign(fleft-fright);


   % V(j,:)=V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));%%�ٶȸ��¹�ʽ
     V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbestgbso-pop(j,:));%%�ٶȸ��¹�ʽ
      


      V(j,find(V(j,:)>Vmax))=Vmax;%�ٶ�Լ��
      V(j,find(V(j,:)<Vmin))=Vmin;

      % ��Ⱥ����
      pop(j,:)=pop(j,:)+(1-k)*V(j,:)+k*Y(j,:);

%�����Ը�
 k=k+0.001;
    if k>1
        k=1;
    end

      pop(j,find(pop(j,:)>popmax))=popmax;%λ��Լ��
      pop(j,find(pop(j,:)<popmin))=popmin;

      %��Ӧ��ֵ����
      fitness(j)=funbas(pop(j,1),pop(j,2));
  end

  for j=1:sizepop
     %��������ֵ����
     if fitness(j)>fitnesspbest(j)
          pbest(j,:)=pop(j,:);
          fitnesspbest(j)=fitness(j);
     end
  %Ⱥ������ֵ����
  if fitness(j)>fitnessgbestgbso
         % gbest=pop(j,:);
      gbestgbso=pop(j,:);
      fitnessgbestgbso=fitness(j);
     
  end
     
  end
gbso(i)=fitnessgbestgbso;
end

%% DE




size=50;%Ⱥ�����
Codel=2;%����ı�������
 MinX(1)=-10;%δ֪����Χ
 MinX(2)=-10;
 MaxX(1)=10;
 MaxX(2)=10;
 G=100;%��������
 F=1.2;%��������[0 2]
 cr=0.8;%��������[0.6 0.9]
 %��ʼ����Ⱥ
 for i=1:1:Codel
 P(:,i)=MinX(i)+(MaxX(i)-MinX(i))*rand(size,1);    
 end
 
 Best=P(1,:);%ȫ�����Ÿ��� ֮�󲻶ϸ���
 
  for i=2:size
     if(fun(P(i,1),P(i,2))>fun(Best(1),Best(2)))
         Best=P(i,:);
     end
  end
  
  
  fi=fun(Best(1),Best(2));%����C���� һ��Ҫ�ǵø���ʼ������������ܷ�
  
  
  %%����ѭ��ֱ�����㾫��Ҫ����ߵ��������ﵽ
   
  for Kg=1:1:G
     time(Kg)=Kg;
    %% �ڶ��� ����
      for i=1:size
          r1=1;r2=1;r3=1;r4=1;%ʹ�ø��������������
          while(r1==r2||r1==r3||r1==r4||r2==r3||r2==r4||r3==r4||r1==i||r2==i||r3==i||r4==i)
            %ҪѰ�ҵ�i�ε����У����ѡ�����������r1 r2  r3��i������Ȳſ���
            r1=ceil(size*rand(1));%��Сƥ�� 
            r2=ceil(size*rand(1));
            r3=ceil(size*rand(1));
            r4=ceil(size*rand(1));
          end
          h(i,:)=P(r1,:)+F*(P(r2,:)-P(r3,:));%p(r2)-p(r3)Ϊ�������
          %h(i,:)=Best+F*(P(r2,:)-P(r3,:));
          for j=1:Codel %����Ƿ�Խ��
              if(h(i,j)<MinX(j))
                  h(i,j)=MinX(j);
              elseif(h(i,j)>MaxX(j)) 
                  h(i,j)=MaxX(j);
              end
          end
          %% ����
        for j=1:Codel
        temper=rand(1);
        if(temper<cr)
            v(i,j)=h(i,j);
        else
            v(i,j)=P(i,j);
        end
        end
        %ѡ��
        if(fun(v(i,1),v(i,2))<fun(P(i,1),P(i,2)))
            P(i,:)=v(i,:);
        end
        if(fun(P(i,1),P(i,2))<fi)
            fi=fun(P(i,1),P(i,2));
            Best=P(i,:);
        end
      end
      Best_f(Kg)=fun(P(i,1),P(i,2));
     % yy=Best_f(Kg);
  end
%% SSA
maxgen=100;   % ��������  
sizepop=50;   %��Ⱥ��ģ
popmax=10;   %λ�ñ߽�
popmin=-10;
k=2;
P_percent = 0.2;
pNum = round( sizepop *  P_percent ); 

for i=1:sizepop
   pop(i,:)=10*rands(1,2);
   
   fitness(i)=fun(pop(i,1),pop(i,2));
    
end

pFit = fitness;                      
[ fMax, bestI ] = max( fitness );      % fMin denotes the global optimum fitness value
bestX = pop( bestI, : );             % bestX denotes the global optimum position corresponding to fMin


%% ��ȸ�����㷨

for t = 1 : maxgen 
 [ ans, sortIndex ] = sort( pFit );% Sort.
     
  [fmax,B]=min( pFit );
   worse= pop(B,:);  
         
   r2=rand(1);
   
   if(r2<0.8)
 
    for i = 1 : pNum                                                   % Equation (3)
         r1=rand(1);
        pop( sortIndex( i ), : ) = pop( sortIndex( i ), : )*exp(-(i)/(r1*maxgen));
        
         pop(sortIndex( i ),find(pop(sortIndex( i ),:)>popmax))=popmax;%λ��Լ��
      pop(sortIndex( i ),find(pop(sortIndex( i ),:)<popmin))=popmin;
      
        fitness(sortIndex( i ))=fun(pop(sortIndex( i ),1),pop(sortIndex( i ),2)); 
    end
  else
  for i = 1 : pNum   
          
  pop( sortIndex( i ), : ) = pop( sortIndex( i ), : )+randn(1)*ones(1,k);
  
    pop(sortIndex( i ),find(pop(sortIndex( i ),:)>popmax))=popmax;%λ��Լ��
      pop(sortIndex( i ),find(pop(sortIndex( i ),:)<popmin))=popmin;
      
  fitness(sortIndex( i ))=fun(pop(sortIndex( i ),1),pop(sortIndex( i ),2)); 
       
  end
     
  pop(sortIndex( i ),find(pop(sortIndex( i ),:)>popmax))=popmax;%λ��Լ��
      pop(sortIndex( i ),find(pop(sortIndex( i ),:)<popmin))=popmin;
  
  
   end
   
   
[ fMMax, bestII ] = min( fitness );      
  bestXX = pop( bestII, : );  
   
   for i = ( pNum + 1 ) : sizepop                     % Equation (4)
     
         A=floor(rand(1,2)*2)*2-1;
         
          if( i>(sizepop/2))
           pop( sortIndex(i ), : )=randn(1)*exp((worse-pop( sortIndex( i ), : ))/(i)^2);
          else
          pop( sortIndex( i ), : )=bestXX+(abs(( pop( sortIndex( i ), : )-bestXX)))*(A'*(A*A')^(-1))*ones(1,k);  

         end  
             
  pop(sortIndex( i ),find(pop(sortIndex( i ),:)>popmax))=popmax;%λ��Լ��
      pop(sortIndex( i ),find(pop(sortIndex( i ),:)<popmin))=popmin;
      
        fitness(sortIndex( i ))=fun(pop(sortIndex( i ),1),pop(sortIndex( i ),2));
        
   end

      c=randperm(numel(sortIndex));
      b=sortIndex(c(1:3));
    for j =  1  : length(b)      % Equation (5)

    if( pFit( sortIndex( b(j) ) )<(fMax) )

        pop( sortIndex( b(j) ), : )=bestX+(randn(1,2)).*(abs(( pop( sortIndex( b(j) ), : ) -bestX)));

        else

        pop( sortIndex( b(j) ), : ) =pop( sortIndex( b(j) ), : )+(2*rand(1)-1)*(abs(pop( sortIndex( b(j) ), : )-worse))/ ( pFit( sortIndex( b(j) ) )-fmax+1e-50);

    end
         
  pop(sortIndex( i ),find(pop(sortIndex( i ),:)>popmax))=popmax;%λ��Լ��
      pop(sortIndex( i ),find(pop(sortIndex( i ),:)<popmin))=popmin;
      
       fitness(sortIndex( i ))=fun(pop(sortIndex( i ),1),pop(sortIndex( i ),2));
    end
   
      for i = 1 : sizepop 
        if ( fitness( i )<pFit( i ) )
            pFit( i ) = fitness( i );
             pop(i,:) = pop(i,:);
        end
        
        if( pFit( i )<fMax )
           fMax= pFit( i );
            bestX =pop( i, : );
         
            
        end
    end
 
 

    ssa(t)=fMax;    
    
    
end


%% �������㷨

N = 50; %��Number of nests(The scale of solution)
D = 10 ; %  Dimensionality of solution
T = 100 ; % Number of iterations
Xmax = 10;
Xmin = -10 ;
Pa = 0.25 ; % Probability of building a new nest(After host bird find exotic bird eggs)
nestPop = rand(N,D)*(Xmax-Xmin)+Xmin ;  % Random initial solutions
i=1;
for t=1:T 
    levy_nestPop =  func_levy(nestPop,Xmax,Xmin) ; % ͨ��Levy flights�����µĽ������
    nestPop = func_bestNestPop(nestPop,levy_nestPop);  %     ���¾ɳ���ѡ��һ����õĳ�
    rand_nestPop = func_newBuildNest(nestPop,Pa,Xmax,Xmin); %������Pa������ĳ�Ѩ��ͨ����ƫ��������ߣ������³�Ѩ
    nestPop = func_bestNestPop(nestPop,rand_nestPop) ; %���¾ɳ���ѡ��һ����õĳ�
    [~,index] = max(func_fitness(nestPop)) ; % Best nests
    trace(t) = func_objValue(nestPop(index,1),nestPop(index,1)) ; 
    diebgb(1)=trace(1);
    
    if diebgb(i)>trace(t)
        diebgb(i+1)=trace(t);
        i=i+1;
    end
    if diebgb(i)<=trace(t)
        diebgb(i+1)=diebgb(i);
        i=i+1;
    end
end
i=i-1;
bgb=min(trace);


%% ����������ͼ
disp('pso')
disp(gbest)
disp(fitnessgbest)

%%%BAS�ָ�������
fbest=-fbest;
fbest_store=-fbest_store;

%%%BSO�ָ�������
fitnessgbestbso=-fitnessgbestbso;
bso=-bso;

%%%gbso�ָ�������
fitnessgbestgbso=-fitnessgbestgbso;
gbso=-gbso;


disp('BAS')
disp(xbest(1));
disp(xbest(2));
disp(fbest);

disp('BSO')
[gbestbso fitnessgbestbso ]

disp('GBSO')
[gbestgbso fitnessgbestgbso ]

disp('DE');
disp(Best);
disp(Best_f(Kg));

disp('SSA');
disp(bestX);
disp(fMax);

disp('bgb');
disp(nestPop(1,1));
disp(nestPop(1,2));
disp(bgb);

% 
% plot3(gbest(1),gbest(2),fitnessgbest,'bo','linewidth',1.5);
% hold on
% plot3(xbest(1),xbest(2),fbest,'r*','linewidth',1.5);
% hold on
% plot3(gbestbso(1),gbestbso(2),fitnessgbestbso,'go','linewidth',1.5);
% hold on
% plot3(gbestgbso(1),gbestgbso(2),fitnessgbestgbso,'ko','linewidth',1.5);
% hold on
% plot3(Best(1),Best(2),Best_f(Kg),'yo','linewidth',1.5);
% hold on
% plot3(bestX(1),bestX(2),fMax,'mo','linewidth',1.5);
% hold on
% plot3(nestPop(1,1),nestPop(1,2),bgb,'co','linewidth',1.5);
figure
plot(yy,'-+r')
hold on
plot(fbest_store,'-*b')
hold on
plot(bso,'-xg')
hold on
plot(gbso,'k')
hold on
plot(Best_f,'-py')
hold on
plot(ssa,'-dm')
hold on
plot(diebgb,'-<c')
legend('PSO','BAS','BSO','RGP-BSO','DE','SSA','CS')

title('���Ÿ�����Ӧֵ','fontsize',12);
xlabel('��������','fontsize',12);
ylabel('��Ӧ��','fontsize',12);





