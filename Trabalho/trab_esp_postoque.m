    dyLdt1 = yL_0(5);
    dyLdt2 = yL_0(6);
    dyLdt3 = yL_0(7); 
    dyLdt4 = yL_0(8);
    dyLdt5 = -(((c_r + c_t)*yL_0(5))/m) + (c_t*yL_0(6))/m + (k_t*(-yL_0(1) + yL_0(2)) + k_r*(-yL_0(1) + y_ext) + c_r*yponto_ext)/m;
    dyLdt6 = (c_t*J_oz*yL_0(5))/(M*(-(M*(D_go*D_go)) + J_oz)) + ((-(M*c_tf*D_fo*D_go) + (c_t + c_tf)*J_oz)*yL_0(6))/(M*(M*(D_go*D_go) - J_oz)) + (c_tf*D_fo*(-(M*D_fo*D_go) + J_oz)*yL_0(7))/(M*(M*(D_go*D_go) - J_oz)) + (M*(D_go*D_go)*(-2*g*M + (2*g*M - S*rho*C_L*((yL_0(10) + u_v)*(yL_0(10) + u_v)))*u_rol*yL_0(3)) - J_oz*(S*rho*C_L*((yL_0(10) + u_v)*(yL_0(10) + u_v)) - 2*(k_t*(-yL_0(1) + yL_0(2)) + k_tf*(yL_0(2) + D_fo*yL_0(3) - yL_0(4)))) + M*D_go*(S*rho*D_po*((yL_0(10) + u_v)*(yL_0(10) + u_v))*(C_L + C_D*yL_0(3)) - 2*D_fo*k_tf*(yL_0(2) + D_fo*yL_0(3) - yL_0(4))))/(2.*M*(M*(D_go*D_go) - J_oz)) + (c_tf*(M*D_fo*D_go - J_oz)*yL_0(8))/(M*(M*(D_go*D_go) - J_oz));
    dyLdt7 = (c_t*D_go*yL_0(5))/(M*(D_go*D_go) - J_oz) + ((c_tf*D_fo - (c_t + c_tf)*D_go)*yL_0(6))/(M*(D_go*D_go) - J_oz) + (c_tf*D_fo*(D_fo - D_go)*yL_0(7))/(M*(D_go*D_go) - J_oz) + (-(S*rho*D_po*((yL_0(10) + u_v)*(yL_0(10) + u_v))*(C_L + C_D*yL_0(3))) + D_go*(S*rho*C_L*((yL_0(10) + u_v)*(yL_0(10) + u_v))*(1 + u_rol*yL_0(3)) - 2*(-(g*M) + k_t*(-yL_0(1) + yL_0(2)) + D_fo*k_tf*yL_0(3) + g*M*u_rol*yL_0(3) + k_tf*(yL_0(2) - yL_0(4)))) + 2*D_fo*k_tf*(yL_0(2) + D_fo*yL_0(3) - yL_0(4)))/(2*M*(D_go*D_go) - 2*J_oz) + (c_tf*(D_fo - D_go)*yL_0(8))/(-(M*(D_go*D_go)) + J_oz);
    dyLdt8 = (c_tf*yL_0(6))/m_f + (c_tf*D_fo*yL_0(7))/m_f - ((c_rf + c_tf)*yL_0(8))/m_f + (k_tf*(yL_0(2) + D_fo*yL_0(3) - yL_0(4)) + k_rf*(-yL_0(4) + y_ext) + c_rf*yponto_ext)/m_f;

    
    
    
    
    
