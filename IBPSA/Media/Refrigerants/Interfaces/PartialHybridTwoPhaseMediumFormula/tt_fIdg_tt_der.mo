within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function tt_fIdg_tt_der
  "Calculates time derivative of tt_fIdg_tt"
  input Real tau "Reduced temperature";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_tt_fIdg_tt "Time derivative of tt_fIdg_tt";

algorithm
  der_tt_fIdg_tt := der_tau*(2*tt_fIdg_tt(tau=tau)+ttt_fIdg_ttt(tau=tau))/tau;

  annotation(Inline=false,
             LateInline=true);
end tt_fIdg_tt_der;
