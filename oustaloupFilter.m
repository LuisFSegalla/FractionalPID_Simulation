function G = oustaloupLuis(frac,N,w_B,w_H)
  pkg load control;
  s = tf('s');
  
  w_B = w_B;
  w_h = w_H;
  
  K = w_H^frac;
  
  num = 1;
  den = 1;
  
  for k = -N:N
    
    wk1 = w_B*( (w_H/w_B)^(((k+N)+(0.5*(1-frac)))/(2*N + 1)) );
    wk2 = w_B*( (w_H/w_B)^(((k+N)+(0.5*(1+frac)))/(2*N + 1)) );
    
    num = num*(s + wk1);
    den = den*(s + wk2);
  endfor
  
  G = ((K*num/den));
end