%--------------      on fait le grand menage      -------------------------
clc
close all
clear all
format short e;
format compact;
dbstop if error;
%--------------------------------------------------------------------------
%-------    Resolution du probleme - u"(x) + u(x) = f(x)   ----------------
%--------------------------------------------------------------------------
                                             % le nombre des noeuds desires
nddl=input('donner le nombre de degre de liberte nddl : '); 
                                                    % le second membre f(x)

f=@(x)(2+x.*(1-x));

                                                  % la solution exacte u(x)

sln = inline('x.*(1-x)');


%--------

                                      % discretisation de l'intervalle [0,1]
xx=linspace(0,1,nddl);  
                                          % definir le pas de chaque maille
for i=1:nddl-1
    h(i)=xx(i+1)-xx(i);
end
h;
                                               % les conditions aux limites
x0=0; 
xL=1;
ua=0; 
ub=0;
                                                   % la matrice de régidité 
A=zeros(nddl-2);

   for i=1:nddl-2
       for j=1:nddl-2
           if(abs(i-j)>1)
               A(i,j)=0;
           elseif(abs(i-j)==0)
               A(i,j)=(1./h(i)) + (1./h(i+1)) + (1./3)*(h(i)+h(i+1));
           else
               A(i,j)=(-1./h(i)) + h(i)./6;
           end
       end
   end

    
                                      % le second membre du probleme primal
                                              
B1=0.5*(h(1)+h(2))*f(xx(2)) - ua*h(1)/4 + ua*(1/h(1));
B_nddl_2=0.5*(h(nddl-2)+h(nddl-1))*f(xx(nddl-1)) - ub*h(nddl-1)/4 + ub*(1/h(nddl-1));

Bb=[];
for i=2:nddl-3
    Bb=[Bb 0.5*(h(i)+h(i+1))*f(xx(i+1))];
    B=[B1 Bb B_nddl_2]; 
end
                                                    % la solution approchee
u=[ua ; A\B' ; ub];
                    % comparaison entre la sol. exacte et la sol. approchee
figure(140)
plot(0:0.01:1,sln(0:0.01:1),'-.K','LineWidth',4)
hold on
plot(xx,u,'c','LineWidth',3)
legend('sol exact','sol approchée pour n=3')
%hold on
title('Graphe de la solution exacte et approchée de u"(x)+u(x)=f(x) avec u(x)=x(1-x)')
pause();
