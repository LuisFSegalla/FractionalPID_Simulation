function [Gmf_Int Gmf_Frac] = compara(xInt,xFrac,Gp,Gf,N,wb,wh,t,vecSim)
    %carrega os pacotes necessários
    pkg load control;
    s = tf('s');
    
    %controlador inteiro
    GcInt = xInt(1) + (xInt(2)/s) + xInt(3)*s;
    
    %controlador fracionário
    f1     = oustaloupLuis(-xFrac(4),N,wb,wh);
    f2     = oustaloupLuis( xFrac(5),N,wb,wh);
    GcFrac = minreal(xFrac(1) + xFrac(2)*f1 + xFrac(3)*f2);
    
    %malha fechada
    Gmf_Int  = minreal(Gp*GcInt  / minreal(1 + Gp*GcInt*Gf));
    Gmf_Frac = minreal(Gp*GcFrac / minreal(1 + Gp*GcFrac*Gf));
    
    %resposta ao degrau
    yInt  = lsim(Gmf_Int,vecSim,t);
    yFrac = lsim(Gmf_Frac,vecSim,t);
    
    %plot
    plot(t,yInt,'linewidth',3,'color','b',t,yFrac,'linewidth',3,'color','r')
    title('Comparação com controlador inteiro e fracionário');
    legend('Ordem Inteira','Ordem Fracionária');
    xlabel('Entrada')
    ylabel('Saída')
    grid on;
       
end