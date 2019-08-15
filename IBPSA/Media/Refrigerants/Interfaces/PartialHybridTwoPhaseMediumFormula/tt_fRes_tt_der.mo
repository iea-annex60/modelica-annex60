within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function tt_fRes_tt_der
  "Calculates time derivative of tt_fRes_tt"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_tt_fRes_tt "Time derivative of tt_fRes_tt";

algorithm
  der_tt_fRes_tt := der_tau*(2*tt_fRes_tt(delta=delta,tau=tau)+
                    ttt_fRes_ttt(delta=delta,tau=tau))/tau +
                    der_delta*ttd_fRes_ttd(delta=delta,tau=tau)/delta;

  annotation(Inline=false,
             LateInline=true);
end tt_fRes_tt_der;
