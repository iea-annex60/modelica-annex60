within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dd_fRes_dd_der
  "Calculates time derivative of dd_fRes_dd"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_dd_fRes_dd "Time derivative of dd_fRes_dd";

algorithm
  der_dd_fRes_dd := der_tau*tdd_fRes_tdd(delta=delta,tau=tau)/tau +
                    der_delta*(2*dd_fRes_dd(delta=delta,tau=tau)+
                    ddd_fRes_ddd(delta=delta,tau=tau))/delta;

  annotation(Inline=false,
             LateInline=true);
end dd_fRes_dd_der;
