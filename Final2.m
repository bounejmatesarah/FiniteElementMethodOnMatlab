%--------------      on fait le grand menage      -------------------------
clc
close all
clear all
format short e;
format compact;
dbstop if error;
%--------------------------------------------------------------------------
%-------     Resolution du probleme - u"(x) + u(x) = f(x)   ---------------
%------------------------     u(x)=x*(x-1)   ------------------------------
%--------------------------------------------------------------------------


                                             % le nombre des noeuds d�sir�s
nddl=input('donner le nombre de degre de liberte nddl : '); 


                     % le second membre f(x)
f=@(x)(-2+(x.*(x-1)));

                                                  % la solution exacte u(x)
sln = inline('x.*(x-1)');



                                     % discretisation de l'intervalle [0,1]
xx=linspace(0,1,nddl);  

s=0;
for i=1:nddl-2
    g=@(x)(abs(x.*(x-1)-(u(i).*(x-xx(i+1))/(xx(i)-xx(i+1)))));
    s=s+quad(g,xx(i),xx(i+1));
end


   % definir le pas de chaque maille cas g�n�ral d'un maillage non uniforme                                       
for i=1:nddl-1
    h(i)=xx(i+1)-xx(i);
end
h;  




                                              % les conditions aux limites
x0=0; 
xL=1;
ua=0; 
ub=0;




                                                   % la matrice de rigidit� 
A=zeros(nddl-2);

  for i=1:nddl-2
       for j=1:nddl-2
           if(abs(i-j)>1)
               A(i,j)=0;
           elseif(abs(i-j)==0)
               A(i,j)=2*(1./h(i)) ;
           else
               A(i,j)=(-1./h(i)) ;
           end
       end
   end

    
   
   
                                      % le second membre du probleme primal
Bb=[];
for i=1:nddl-2
    Bb=[Bb 0.5*(h(i)+h(i+1))*f(xx(i+1))];
    B=[Bb ]; 
end



                                                    % la solution approchee
u=[ua ; A\B' ; ub];



                    % comparaison entre la sol.exacte et la sol. approchee
figure(140)
plot(0:0.01:1,sln(0:0.01:1),'-.K','LineWidth',3)
hold on
plot(xx,u,'o-','Color',[ 8.7059e-001  4.9020e-001            0],'LineWidth',2)
legend('\bf sol exacte','\bf sol approch�e pour n=3')
%hold on
xlabel('\fontname{Times}\fontsize{14} \bf x','Color',[ 8.7059e-001  4.9020e-001            0])
ylabel('\fontname{Times}\fontsize{14} \bf u(x)','Color',[ 8.7059e-001  4.9020e-001            0]);
grid

set(gcf,'Color',[ 9.5294e-001  8.7059e-001  7.3333e-001])
whitebg([ 9.5294e-001  8.7059e-001  7.3333e-001]);
