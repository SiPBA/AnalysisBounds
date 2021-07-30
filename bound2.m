% Implementation of (non)-linear upper bounds
% Please cite
% https://www.sciencedirect.com/science/article/abs/pii/S0031320319301402
% J.M.Gorriz,  et al. 
% On the computation of distribution-free performance bounds: Application to small sample sizes in neuroimaging. 
% Pattern Recognition 93, 1-13, 2019.

clear;
clc;
%l=120; % Sample Size of your experiment (number of training features)
l=500;
eta=0.05; % level of signiifcance (...with probability at least...)
r=1:1:10; % Complecity of the classifier 
dim=10;   % Maximun number of predictors or feature dimensions 
D=zeros(numel(r),dim);
for i=1:numel(r)
    D(i,max(r(i))+1:end)=max(r(i))+1:1:dim; %feature dimension
end
%d=3;

BoundGr=zeros(numel(r),dim);
Cldr=zeros(numel(r),dim);

for s=1:numel(r)
    
    for j=1:dim
        if D(s,j)>r(s)
            jnew=nchoosek(D(s,j),r(s));
            for k=1:min(jnew,l-1)
                Cldr(s,j)=Cldr(s,j)+nchoosek(l-1,k-1);
            end
            Cldr(s,j)= 2*Cldr(s,j);
            BoundGr(s,j)=sqrt(log(Cldr(s,j)/eta)/(2*l));
        end
    end
    
        
end

%2D
figure1=figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(D',BoundGr','*-')
clear leyenda
for i=1:numel(r)
    leyenda{i}=['complexity ' int2str(r(i))];
end
legend(leyenda)
legend1 = legend(axes1,'show');

% Create ylabel
ylabel('\Delta_n(Z^n)','FontWeight','bold');

% Create xlabel
xlabel('dimension','FontWeight','bold');

box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontWeight','bold');

set(legend1,...
    'Position',[0.670535714285714 0.554761890874429 0.124107142857143 0.282142857142857]);

