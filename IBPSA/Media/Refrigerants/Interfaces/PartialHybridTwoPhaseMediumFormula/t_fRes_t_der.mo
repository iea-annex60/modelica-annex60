within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function t_fRes_t_der
  "Calculates time derivative of t_fRes_t"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_t_fRes_t "Time derivative of t_fRes_t";

algorithm
  der_t_fRes_t := der_tau*(t_fRes_t(delta=delta,tau=tau) +
                  tt_fRes_tt(delta=delta,tau=tau))/tau + der_delta*
                  td_fRes_td(delta=delta,tau=tau)/delta;

  annotation(Inline=false,
             LateInline=true);
end t_fRes_t_der;
