clc; clear; close all;
pkg load control;
pkg load statistics;
s = tf('s');

%---------------------------------------------------------%
%Parâmetros do Filtro
N  = 3;
wb = 10^-3;
wh = 10^3;
t = 0:0.2:10;

[num,den] = padecoef(1,5);
delay = tf(num,den);


M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;
q = (M+m)*(I+m*l^2)-(m*l)^2;
s = tf('s');

P_cart = (((I+m*l^2)/q)*s^2 - (m*g*l/q))/(s^4 + (b*(I + m*l^2))*s^3/q - ((M + m)*m*g*l)*s^2/q - b*m*g*l*s/q);

P_pend = (m*l*s/q)/(s^3 + (b*(I + m*l^2))*s^2/q - ((M + m)*m*g*l)*s/q - b*m*g*l/q);

%planta
G_p = (-s+1)/ ((s^2+1)*(s+1))
G_f = 1

%Com realimentação não unitária
costFunction = @(x,G_p) ITAEOustaloup_Control(x,G_p,G_f,N,wb,wh);
%costFunction = @(x,G_p) ITAE_Control(x,G_p,G_f);


nVar         = 5;                %número de variáveis de decisão
varSize      = [1 nVar];         %Matriz de variáveis de decisÃ£o
varMin       = [0.01 0.01 0.01 0.01 0.01];%mínimo valor que uma variável pode assumir
varMax       = [  10   10   10    2    2];%máximo valor que uma variável pode assumir
maxVelocity  = 0.1*(varMax-varMin);
minVelocity  = -maxVelocity;

%Constriction coeficients | Quando usamos este método não é necessário utilizar o coeficiente de amortecimento para w
kappa = 1;
phi1  = 2.05;
phi2  = 2.05;
phi   = phi1 + phi2;
chi   = 2*kappa / abs(2-phi-sqrt(phi*phi -4*phi));


%%Parâmetros do problema
maxIteration = 200;%máximo de iterações
nPop         = 25;%tamanho da população
w            = chi;%coeficiente de inércia
wDamp        = 1;%Fator de amortecimento do coeficiente de inércia
c1           = chi*phi1;%coeficiente pessoal de aceleração
c2           = chi*phi2;%coeficiente social de aceleração

%Inicialização

%Template da struct particle
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost     = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost     = [];

%cria o vetor da população
particle = repmat(empty_particle,nPop,1);%vetor de empty_particle

%Inicializa o melhor resultado global
GlobalBest.Cost     = inf;


%Inicializa a população
for i=1:nPop
  %gera uma solução aleatória dentro dos limites do domínio
  particle(i).Position = unifrnd(varMin,varMax,varSize);
  
  %inicializa a velocidade 
  particle(i).Velocity = zeros(varSize);
  
  %Avalia a solução encontrada
  particle(i).Cost = costFunction(particle(i).Position,G_p); %ITAE
  
  %atualiza a melhor resposta particular
  particle(i).Best.Position = particle(i).Position;
  particle(i).Best.Cost     = particle(i).Cost;
  
  %atualiza a melhor resposta global
  if(particle(i).Best.Cost < GlobalBest.Cost)
    GlobalBest = particle(i).Best;
  endif
  
endfor

bestCosts = zeros(maxIteration,1);%vetor para guardar oos melhores valores em cada iterações

tic%inicia a medição de tempo do programa

%%Iterações (main loop)
for it=1:maxIteration
    
  for i=1:nPop
    %atualiza a velocidade
    particle(i).Velocity = w*particle(i).Velocity + c1*rand(varSize).*(particle(i).Best.Position - particle(i).Position) + c2*rand(varSize).*(GlobalBest.Position - particle(i).Position);
    
    %avalia a velocidade dentro dos limites do problema
    particle(i).Velocity = max(particle(i).Velocity,minVelocity);
    particle(i).Velocity = min(particle(i).Velocity,maxVelocity);
    
    %atualiza a posição
    particle(i).Position = particle(i).Position + particle(i).Velocity;
    
    %avalia a posição dentro dos limites do problema
    particle(i).Position = max(particle(i).Position,varMin);
    particle(i).Position = min(particle(i).Position,varMax);
    
    
    %avalia a nova posição
    particle(i).Cost = costFunction(particle(i).Position,G_p); %ITAE
    
    %atualiza o melhor resultado pessoal
    if(particle(i).Cost < particle(i).Best.Cost)
      particle(i).Best.Position = particle(i).Position;
      particle(i).Best.Cost     = particle(i).Cost;
      
      %atualiza o melhor resultado global
      if(particle(i).Best.Cost < GlobalBest.Cost)
        GlobalBest = particle(i).Best;
      endif
      
    endif
    
  endfor
  %atualiza o melhor valor de cada iteração
  bestCosts(it) = GlobalBest.Cost;
  disp(['Iteração: ' num2str(it), ' Melhor custo = ' num2str(bestCosts(it))])
  
  %amortecendo o coeficiente de inércia
  w = w * wDamp;
  
endfor

toc%finaliza a medição de tempo do programa

%%Resultados
figure;
semilogy(bestCosts, 'linewidth', 3);
xlabel('Iterações');
ylabel('Custo');
title('Função de custo com o passar das iterações');
grid on;