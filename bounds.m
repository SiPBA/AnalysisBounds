% Implementation of linear upper-bounds 
% Please cite
% https://www.sciencedirect.com/science/article/abs/pii/S0031320319301402
% J.M.Gorriz,  et al. 
% On the computation of distribution-free performance bounds: Application to small sample sizes in neuroimaging. 
% Pattern Recognition 93, 1-13, 2019.

clear;
close all;
clc;

l=[250:1000:41486]; %sample size
eta=0.05; % level of signiifcance (...with probability at least...)
d=1:9:100; % Number of predictors or feature dimensions 
h=d+1; % VC dimension

%% Vapnik´ Bound

BoundVC=zeros(numel(l),numel(d));
for i=1:numel(d)
    for j=1:numel(l)
    BoundVC(j,i)=sqrt(abs((h(i).*(log(2*l(j)./h(i))+1)-log(eta/4))./l(j)));
    end
end

% Bound number 1
BoundG=zeros(numel(l),numel(d));
Cld=zeros(numel(l),numel(d));

for s=1:numel(l)
    
    for j=1:numel(d)
        for k=1:j
            Cld(s,j)=Cld(s,j)+nchoosek(l(s)-1,k-1);
        end
        Cld(s,j)= 2*Cld(s,j);
        BoundG(s,j)=sqrt(log(Cld(s,j)/eta)/(2*l(s)));
    end
    
        
end

% %Calculo de Phi(n,d)
% BoundPhi=zeros(numel(l),numel(d));
% CldPhi=zeros(numel(l),numel(d));
% 
% for s=1:numel(l)
%     
%     for j=1:numel(d)
%         for k=1:j
%             CldPhi(s,j)=CldPhi(s,j)+nchoosek(l(s),k-1);
%         end
%         %Cld(s,j)= 2*Cld(s,j);
%         BoundPhi(s,j)=sqrt(log(CldPhi(s,j)/eta)/(2*l(s)));
%     end
%     
%         
% end

%% Bound Number 2 
%number of functions with k zeroes
BoundGZ=zeros(numel(l),numel(d));
Cld=zeros(numel(l),numel(d));
K=d-1;

for s=1:numel(l)
    
    for j=1:numel(d)
        for z=K
        for k=1:j-z
            Cld(s,j)=Cld(s,j)+nchoosek(l(s),z)*nchoosek(l(s)-1-z,k-1);
        end
        end
        Cld(s,j)= 2*Cld(s,j);
        BoundGZ(s,j)=sqrt(log(Cld(s,j)/eta)/(2*l(s)));
    end
    
        
end



%%%%%%%%%%%%%%%%%PLOTS and GRAPHS%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D
figure
[X,Y]=meshgrid(l,d);
ax1=mesh(X,Y,BoundVC','FaceColor','b');
hold on
ax2=mesh(X,Y,BoundG','FaceColor','g');
ax3=mesh(X,Y,BoundGZ','FaceColor','c');
%ax4=mesh(X,Y,BoundPhi','FaceColor','y');
xlabel('Sample Size l')
ylabel('Dimension d')
zlabel('Upper bound')
legend('Vapnik´s bound \gamma_{VC}','Proposed Bound \gamma_{emp} with N(l,d)','Proposed Bound with \gamma_{emp} with Q(l,d)')
hold off

%2D
figure1=figure;
hold on;


plot(l,BoundG)
%clear leyenda
for i=1:numel(d)
    leyenda{i}=['dim ' int2str(d(i))];
end
legend(leyenda)
%legend1 = legend(axes1,'show');

% Create ylabel
ylabel('\Delta_n(Z^n)','FontWeight','bold');

% Create xlabel
xlabel('samples','FontWeight','bold');

grid on;

% Create arrow
annotation(figure1,'arrow',[0.221626984126984 0.329861111111111],...
    [0.185762778505898 0.603211009174312],'Color',[0 0 1],'LineWidth',3,...
    'LineStyle','--',...
    'HeadStyle','cback2');

title('Upper bound analysis')