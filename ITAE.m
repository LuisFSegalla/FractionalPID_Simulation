function out=ITAE_Control(x,G_ma,Gf)
  pkg load control;
  s = tf('s');
  
  G_c = (x(3)*s*s + x(1)*s + x(2)) / s;
  G_MF = minreal((G_ma*G_c) / (1 + G_ma*G_c*Gf));
  
  t = 0:0.2:10;
  respDegrau = step(G_MF,t);
  degrau     = ones(length(t),1);
  erro       = abs(degrau - respDegrau);
  tt         = t';
  eat        = tt.*erro;
  ITAE       = trapz(tt,eat);  
  out        = ITAE;
end