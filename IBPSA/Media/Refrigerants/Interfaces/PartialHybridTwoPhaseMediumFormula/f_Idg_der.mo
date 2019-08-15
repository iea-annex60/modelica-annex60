within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function f_Idg_der "Calculates time derivative of f_Idg"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced f_Idg";
  output Real der_f_Idg "Time derivative of alpha_0";

algorithm
  der_f_Idg := der_delta/delta + der_tau*t_fIdg_t(tau=tau)/tau;

  annotation(Inline=false,
             LateInline=true);
end f_Idg_der;
