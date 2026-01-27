function model = CAF2_with_denom_v2()
    % CAF2_with_denom_v2 provides the GenSSI implementation of model for cancer dynamics
    % with CAFs 

    % Symbolic variables
	syms alpha2 d3 k2 alpha3 u1 u2 u3
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12

    % Parameters
	model.sym.p = [k2; alpha2; d3; alpha3];

    % State variables
	model.sym.x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12];

    % Control vectors (g)
	model.sym.g = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, - 1*u1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 1*u1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, - 1*u2, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 1*u2, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, - 1*u3, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 1*u3, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    % Autonomous dynamics (f)
	model.sym.xdot = [x1*(1-(0.001*x1))*(3*(x3/(10 + x3)) + k2*(x4/(alpha2 + x4)))
                       5*(1/(1+d3*x4)) - 0.0042*x2 + 0.0417*x3
                       -0.0125*x3 - 0.0417*x3
                       1*x4*(1-(0.001*x4))*(x1/(alpha3 + x1))
                       x5*(1-(0.001*x5))*(3*(x7/(10 + x7)) + k2*(x8/(alpha2 + x8)))
                       5*(1/(1+d3*x8)) - 0.0042*x6 + 0.0417*x7
                       -0.0125*x7 - 0.0417*x7
                       1*x8*(1-(0.001*x8))*(x5/(alpha3 + x5))
                       x9*(1-(0.001*x9))*(3*(x11/(10 + x11)) + k2*(x12/(alpha2 + x12)))
                       5*(1/(1+d3*x12)) - 0.0042*x10 + 0.0417*x11
                       -0.0125*x11- 0.0417*x11
                       1*x12*(1-(0.001*x12))*(x9/(alpha3 + x9))];

    % Initial conditions
	model.sym.x0 = [0.5;0.1;0.01;1.5;0.5;0.1;0.01;1.5;0.5;0.1;0.01;1.5];

    % Observables    
	model.sym.y = [x1;x2;x5;x9];
end
