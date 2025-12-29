function model = CAF2()
    % HIV provides the GenSSI implementation of model for HIC dynamics
    % described by
    % 
    %    Meshkat et al. (2014). On finding and using identifiable parameter
    %    combinations in nonlinear dynamic Systems Biology models and
    %    COMBOS: a novel Web implementation, PLoS ONE, 9, e110261.

    % Symbolic variables
	syms alpha1 betaa k1_hat db x20 x30 u1 u2 u3
    syms x1 x2 x3 x5 x6 x7 x9 x10 x11

    % Parameters
	model.sym.p = [alpha1; betaa; k1_hat; db; x20; x30];

    % State variables
	model.sym.x = [x1 x2 x3 x5 x6 x7 x9 x10 x11];

    % Control vectors (g)
	model.sym.g = [0, 0, 0, 0, 0, 0, 0, 0, 0
        0, - db*u1, 0, 0, 0, 0, 0, 0, 0
        0, db*u1, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, - db*u2, 0, 0, 0, 0
        0, 0, 0, 0, db*u2, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0, 0, - db*u3, 0
        0, 0, 0, 0, 0, 0, 0, db*u3, 0];

    % Autonomous dynamics (f)
	model.sym.xdot = [x1*(1-(0.001*x1))*(k1_hat*(x3/(alpha1 + x3)));
                       betaa - 0.0042*x2 + 0.0417*x3;
                       -0.0125*x3 - 0.0417*x3;
                       x5*(1-(0.001*x5))*(k1_hat*(x7/(alpha1 + x7)));
                       betaa - 0.0042*x6 + 0.0417*x7;
                       -0.0125*x7 - 0.0417*x7;
                       x9*(1-(0.001*x9))*(3*(x11/(alpha1 + x11)) );
                       betaa - 0.0042*x10 + 0.0417*x11;
                       -0.0125*x11 - 0.0417*x11                       ];

    % Initial conditions
	model.sym.x0 = [0.5;x20;x30;0.5;x20;x30;0.5;x20;x30];

    % Observables    
	model.sym.y = [x1;x2;x5;x9];
end
