function out=ITAEOustaloup_Control(x,G_ma,Gf,N,wb,wh)
  pkg load control;
  s = tf('s');
  
  ousta1 = oustaloupLuis(-x(4),N,wb,wh);
  ousta2 = oustaloupLuis(x(5),N,wb,wh);
  
  G_c = minreal(x(1) + x(2)*ousta1 + x(3)*ousta2);
  
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