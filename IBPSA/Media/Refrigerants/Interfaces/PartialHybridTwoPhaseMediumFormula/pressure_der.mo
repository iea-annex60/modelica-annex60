within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function pressure_der
  "Calculates time derivative of pressure calculated by thermodynamic
  state record"
  input ThermodynamicState state "Thermodynamic state";
  input ThermodynamicState der_state
    "Time derivative of thermodynamic state";
  output Real der_p "Time derivative of pressure";

algorithm
  der_p := der_state.p;

  annotation(Inline=true);
end pressure_der;
