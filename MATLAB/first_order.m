function dydt = first_order(t,y,u)
    tau = 5.0;
    K = 2.0;
    
    dydt = (-y+K*u)/tau;
end
    