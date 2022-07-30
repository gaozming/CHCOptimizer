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
%               Cock-hen-chicken Optimizer: A Nature-inspired       %
%               Algorithm for Global Optimization                   %
%               DOI: https://doi.org/10.1155/2019/2981282           %
%                                                                   %
%___________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run ALO: [Best_score,Best_pos,cg_curve]=ALO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________

clear all 
clc

SearchAgents_no=50; % Number of search agents

Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=200; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[~,cg_curve,~,Best_pos,~]=CHC(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Test function')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
grid off

%Draw objective space
subplot(1,2,2);
semilogy(cg_curve,'Color','r')
title('Convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight;
grid off;
box on;
legend('CHC');

display(['The best solution obtained by ALO is : ', num2str(Best_pos(end,:))]);
display(['The best optimal value of the objective funciton found by CHC is : ', num2str(cg_curve(end))]);

        



