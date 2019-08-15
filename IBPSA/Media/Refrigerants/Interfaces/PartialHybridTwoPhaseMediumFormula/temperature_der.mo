within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_der
  "Calculates time derivative of temperature calculated by thermodynamic
  state record"
  input ThermodynamicState state "Thermodynamic state";
  input ThermodynamicState der_state
    "Time derivative of thermodynamic state";
  output Real der_T "Time derivative of temperature";

algorithm
  der_T := der_state.T;

  annotation(Inline=true);
end temperature_der;
