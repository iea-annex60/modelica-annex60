within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function td_fRes_td_der
  "Calculates time derivative of td_fRes_td"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_td_fRes_td "Time derivative of td_fRes_td";

algorithm
  der_td_fRes_td := der_tau*(td_fRes_td(delta=delta,tau=tau)+ttd_fRes_ttd(
                    delta=delta,tau=tau))/tau + der_delta*(td_fRes_td(
                    delta=delta,tau=tau) + tdd_fRes_tdd(
                    delta=delta,tau=tau))/delta;

  annotation(Inline=false,
             LateInline=true);
end td_fRes_td_der;
