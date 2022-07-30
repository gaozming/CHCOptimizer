%___________________________________________________________________%
%  Cock-Hen-Chicken Optimizer                 %
%  source codes version 1.0                                         %
%                                                                   %
%  Developed in MATLAB R2021b(9.11)                                 %
%                                                                   %
%  Author and programmer: Zheng-Ming GAO                            %
%                                                                   %
%         e-Mail: gaozming@jcut.edu.cn                              %
%                 gaozming@sina.com                                 %
%                                                                   %
%       Homepage: http://218.199.48.7:8080/zhuye/201807001          %
%                                                                   %
%   Main paper: Zheng-Ming GAO, Juan ZHAO.                          %
%               An Improved Grey Wolf Optimization Algorithm with   %
%               Variable Weights, Computational Intelligence and    %
%               Neuroscience,2019,2981282.                          %
%               DOI: https://doi.org/10.1155/2019/2981282           %
%                                                                   %
%___________________________________________________________________%

% Cock-Hen-Chicken optimization algorithm
function [fitness_all,Convergence_curve,position_all,position_best,...
    position_first]=CHC(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

if nargin==0
    pathName='unimodalScalableFcn';
    % pathName='multimodalScalableFcn';
    Function_names={
        % unimodal benchmark functions
    %     'Ackley1Fcn';
%         'ExponentialFcn';
%         'PowellSumFcn';
        'SarganFcn';
    %     'Schwefel2p20Fcn';
    %     'Schwefel2p21Fcn';
%         'SphereFcn';
%         'StepFcn';
    %     'Step2Fcn';
    %     'Step3Fcn'
    % };
        % multimodal benchmark functions
    %     'Alpine1Fcn';
    %     'CosineMixtureFcn'
        };
    % adding path
    addpath(genpath(['./',pathName,'/']));
    SearchAgents_no=50;
    Max_iter=100;
    lb=0;
    ub=100;
    dim=10;
    fobj=str2func(Function_names{1});
end

% parameters
cocks=10;
hens=5;
chickens=SearchAgents_no-cocks-hens;

%Initialize the positions of search agents
% Positions=initialization(SearchAgents_no,dim,ub,lb);
positions_cocks=initialization(cocks,dim,ub,lb);
positions_hens=initialization(hens,dim,ub,lb);
positions_chickens=initialization(chickens,dim,ub,lb);
Positions=[positions_cocks;positions_hens;positions_chickens];
% fitness
fit_cocks=zeros(1,cocks);
fit_hens=zeros(1,hens);
fit_chickens=zeros(1,chickens);
fitness=zeros(1,SearchAgents_no);
for i =1:SearchAgents_no
    fitness(i)=fobj(Positions(i,:));
end

Convergence_curve=zeros(1,Max_iter);


it=0;% Loop counter
[best_fit,index]=min(fitness);
best_position=Positions(index,:);
[worst_fit,index]=max(fitness);
worst_position=Positions(index,:);
p1=0.97;

% initializing output
fitness_all=zeros(Max_iter,SearchAgents_no);
Convergence_curve=zeros(1,Max_iter);
position_all=zeros(Max_iter,SearchAgents_no,dim);
position_best=zeros(Max_iter,dim);
position_first=zeros(Max_iter,dim);

% Main loop
while it<Max_iter
    for i=1:SearchAgents_no
        
       % Return back the search agents that go beyond the boundaries of the
       % search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+...
            lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness(i)=fobj(Positions(i,:));
        
        % Update Alpha, Beta, and Delta
        if fitness(i)<best_fit 
            best_fit=fitness(i); % Update alpha
            best_position=Positions(i,:);
        end        
    end

    % constructs the cocks, hens, and chickens
    [fit_sorted,fit_index]=sort(fitness,'ascend');
    position_cocks=Positions(fit_index(1:cocks),:);
    fit_cocks=fit_sorted(1:cocks);
    position_hens=Positions(fit_index(cocks+1:cocks+hens),:);
    fit_hens=fit_sorted(cocks+1:cocks+hens);
    position_chickens=Positions(fit_index(cocks+hens+1:end),:);
    fit_chickens=fit_sorted(cocks+hens+1:end);
    
    a=2-it*2.0/Max_iter; % a decreases linearly fron 2 to 0
    
    % update the cocks
    position_cocks_current=position_cocks;
    for i=1:cocks
        if fit_cocks(i)==best_fit
            for j=1:dim
                % best one , random walk
                position_cocks(i,j)=Positions(i,j)+unifrnd(-1,1)*0.8;
            end
        else
            LF=levyFlight(dim);
            for j =1:dim
                % Levy fligh
                position_cocks(i,j)=best_position(j)*LF(j);
            end
        end
        fit_cocks_new=fobj(position_cocks(i,:));
        % memory saving
        if fit_cocks_new<fit_cocks(i)
            fit_cocks(i)=fit_cocks_new;
        else
            position_cocks(i,:)=position_cocks_current(i,:);
        end
    end

    % update the hens
    position_hens_current=position_hens;
    for i=1:hens
        for j =1:dim
            if rand<0.5
                % towards the best
                position_hens(i,j)=best_position(j)+a*rand*(best_position(j)-position_hens(i,j));
            else
                % towards a random cock
                k=unidrnd(cocks);
                position_hens(i,j)=position_cocks(k,j)+a*rand*(position_cocks(k,j)-position_hens(i,j));
            end
        end
        fit_hens_new=fobj(position_hens(i,:));
        % memory saving
        if fit_hens_new<fit_hens(i)
            fit_hens(i)=fit_hens_new;
        else
            position_hens(i,:)=position_hens_current(i,:);
        end
    end
    % update the chickens
    for i=1:SearchAgents_no-cocks-hens
        for j=1:dim
            if rand<0.5
                if rand<0.5
                    % towards the worst
                    position_chickens(i,j)=position_chickens(i,j)+a*rand*(position_chickens(i,j)-worst_position(j));
                else
                    % towards the best
                    position_chickens(i,j)=position_chickens(i,j)+a*rand*(position_chickens(i,j)-best_position(j));
                end
            else
                if rand<p1
                    % following the hens
                    k=unidrnd(hens);
                    position_chickens(i,j)=position_hens(k,j)+a*rand*(position_hens(k,j)-position_chickens(i,j));
                else
                    % re-initialize
                    if length(ub)==1
                        position_chickens(i,j)=(ub-lb)*rand+lb;
                    else
                        position_chickens(i,j)=(ub(j)-lb(j))*rand+lb(j);
                    end
                end
            end
        end
    end

    % reorganize the swarms
    Positions=[position_cocks;position_hens;position_chickens];
    fitness=[fit_cocks,fit_hens,fit_chickens];
    
    it=it+1;
    % updating the output
    fitness_all(it,:)=fitness(:);
    Convergence_curve(it)=best_fit;
    position_all(it,:,:)=Positions;
    position_best(it,:)=best_position;
    position_first(it,:)=Positions(1,:);

    if mod(it-1,50)==0
        sprintf("CHC optimizer for %d and the best result is %f",it-1,best_fit);
    end
end

%% run alone
if nargin==0
    % removing unimodal benchmark functions
    % moving parallel
    % delete(gcp('nocreate'));
    rmpath(genpath(['./',pathName,'/']));
    toc;
    semilogy(1:Max_iter,Convergence_curve);
end