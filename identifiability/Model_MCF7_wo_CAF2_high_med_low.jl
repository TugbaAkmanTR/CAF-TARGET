using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = x1(t)*(1-(0.001*x1(t)))*(k1_hat*(x3(t)/(alpha1 + x3(t)))),
    x2'(t) = beta - 0.0042*x2(t) - db*u1(t)*x2(t) + 0.0417*x3(t),
    x3'(t) = -0.0125*x3(t) + db*u1(t)*x2(t) - 0.0417*x3(t),
    x4'(t) = x4(t)*(1-(0.001*x4(t)))*(k1_hat*(x6(t)/(alpha1 + x6(t)))),
    x5'(t) = beta - 0.0042*x5(t) - db*u2(t)*x5(t) + 0.0417*x6(t),
    x6'(t) = -0.0125*x6(t) + db*u2(t)*x5(t) - 0.0417*x6(t),
    x7'(t) = x7(t)*(1-(0.001*x7(t)))*(k1_hat*(x9(t)/(alpha1 + x9(t)))),
    x8'(t) = beta - 0.0042*x8(t) - db*u3(t)*x8(t) + 0.0417*x9(t),
    x9'(t) = -0.0125*x9(t) + db*u3(t)*x8(t) - 0.0417*x9(t),
    y1(t) = x1(t),
    y2(t) = x2(t),
)

println(assess_identifiability(ode,known_ic=[x1,x4,x7]))
#println(assess_identifiability(ode,known_ic=[x1,x2,x3,x4,x5,x6,x7,x8,x9]))

# Result
# Info: Assessing identifiability with known initial conditions concluded in 0.5816821 seconds
#OrderedCollections.OrderedDict{Nemo.QQMPolyRingElem, Symbol}(
#x1(0) => :globally, 
#x2(0) => :globally, 
#x3(0) => :globally, 
#x4(0) => :globally, 
#x5(0) => :nonidentifiable, 
#x6(0) => :nonidentifiable, 
#x7(0) => :globally, 
#x8(0) => :nonidentifiable, 
#x9(0) => :nonidentifiable, 
#alpha1 => :globally, 
#beta => :globally, 
#db => :globally, 
#k1_hat => :globally)

