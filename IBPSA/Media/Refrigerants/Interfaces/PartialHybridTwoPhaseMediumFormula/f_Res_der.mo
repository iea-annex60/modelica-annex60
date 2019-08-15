within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function f_Res_der "Calculates time derivative of f_Res"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_f_Res "Time derivative of f_Res";

algorithm
  der_f_Res := der_delta*d_fRes_d(delta=delta,tau=tau)/delta +
               der_tau*t_fRes_t(delta=delta,tau=tau)/tau;

  annotation(Inline=false,
             LateInline=true);
end f_Res_der;
