within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function density_der
  "Calculates time derivative of density calculated by thermodynamic
  state record"
  input ThermodynamicState state "Thermodynamic state";
  input ThermodynamicState der_state
    "Time derivative of thermodynamic state";
  output Real der_d "Time derivative of density";

algorithm
  der_d := der_state.d;

  annotation(Inline=true);
end density_der;
