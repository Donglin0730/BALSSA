    % -----------------------------------------------------------------------------------------------------------
% Sparrow Search algorithm (SSA) (demo)
% Programmed by Jian-kai Xue    
% Updated 21 Mar, 2020.                     
%
% This is a simple demo version only implemented the basic         
% idea of the SSA for solving the unconstrained problem.    
% The details about SSA are illustratred in the following paper.    
% (To cite this article):                                                
%  Jiankai Xue & Bo Shen (2020) A novel swarm intelligence optimization
% approach: sparrow search algorithm, Systems Science & Control Engineering, 8:1, 22-34, DOI:
% 10.1080/21642583.2019.1708830


clear all 
clc
h=1;
i=1;
while(i<=10)
Function_name='F12'; 
M=1000;     % Maximum numbef of iterations
pop=30;  
[c,d,dim,fobj] = Get_Functions_details(Function_name);
 

[fMin,bestX,Convergence_curve] =BALSSA(pop,M,c,d,dim,fobj);

display(['The best optimal value of the objective funciton found by SSA is : ', num2str(fMin)]);
A(i,:)=fMin;
B(i,:)=Convergence_curve;

i=i+1;
end
aa(h,:)=min(A);
a=mean(B);
h=h+1;
aa(h,:)=mean(A);
h=h+1;
aa(h,:)=std(A);%BALSSA
semilogy(a,'color',[0.8500 0.3250 0.0980],'Linewidth',2)
hold on;
legend('BALSSA')
axis tight
grid on
box on

