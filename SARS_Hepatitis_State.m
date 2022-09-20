function dx = SARS_Hepatitis_State(t,y)

   global Lambda mu2 beta_C2 beta_H2 beta_CH2 varthet_12 varthet_22 alpha_12 alpha_22 eta_C2 eta_H2 eta_CH2 zet_12 zet_22 xi_C2 xi_H2 xi_CH2
   
dx = zeros(7,1);

dx(1)= Lambda - beta_C2*(varthet_12*y(2) + y(3))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(1) - beta_H2*(varthet_22*y(4) + y(5))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(1) - (beta_CH2*y(6))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(1) - mu2*y(1); 
dx(2) = beta_C2*(varthet_12*y(2) + y(3))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(1) - (alpha_12 + mu2)*y(2); 
dx(3) = alpha_12*y(2) - (xi_C2 + eta_C2 + mu2)*y(3) - zet_12*beta_H2*(varthet_22*y(4) + y(5))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(3);
dx(4) = beta_H2*(varthet_22*y(4) + y(5))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(1) - (alpha_22 + mu2)*y(4); 
dx(5) = alpha_22*y(4) - (xi_H2 + eta_H2 + mu2)*y(5) - zet_22*beta_C2*(varthet_12*y(2) + y(3))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(5); 
dx(6) = (beta_CH2*y(6))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(1) + zet_12*beta_H2*(varthet_22*y(4) + y(5))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(3) + zet_22*beta_C2*(varthet_12*y(2) + y(3))/(y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7))*y(5) - (xi_CH2 + eta_CH2 + mu2)*y(6); 
dx(7) = xi_C2*y(3) + xi_H2*y(5) + xi_CH2*y(6) - mu2*y(7);