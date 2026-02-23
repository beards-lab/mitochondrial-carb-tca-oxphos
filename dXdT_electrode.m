function f = dXdT_electrode(t,Jel,tsim,Jox)
tau = 5; % seconds
f = (interp1(tsim,Jox,t)- Jel)/tau;