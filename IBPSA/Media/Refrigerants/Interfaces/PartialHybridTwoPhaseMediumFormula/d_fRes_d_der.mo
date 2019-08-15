within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function d_fRes_d_der
  "Calculates time derivative of d_fRes_d"
  input Real delta "Reduced density";
  input Real tau "Reduced temperature";
  input Real der_delta "Time derivative of reduced density";
  input Real der_tau "Time derivative of reduced temperature";
  output Real der_d_fRes_d "Time derivative of d_fRes_d";

algorithm
  der_d_fRes_d := der_tau*td_fRes_td(delta=delta,tau=tau)/tau +
                  der_delta*(d_fRes_d(delta=delta,tau=tau)+
                  dd_fRes_dd(delta=delta,tau=tau))/delta;

  annotation(Inline=false,
             LateInline=true);
end d_fRes_d_der;
