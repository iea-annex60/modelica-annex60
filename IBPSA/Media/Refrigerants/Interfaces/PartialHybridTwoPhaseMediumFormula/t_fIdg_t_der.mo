within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function t_fIdg_t_der
  "Calculates time derivative of t_fIdg_t"
  input Real tau "Reduced temperature";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_t_fIdg_t "Time derivative of t_fIdg_t";

algorithm
  der_t_fIdg_t := der_tau*(t_fIdg_t(tau=tau)+ tt_fIdg_tt(tau=tau))/tau;

  annotation(Inline=false,
             LateInline=true);
end t_fIdg_t_der;
