within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_der
  "Calculates time derivative of specific enthalpy calculated by
  thermodynamic state record"
  input ThermodynamicState state "Thermodynamic state";
  input ThermodynamicState der_state
    "Time derivative of thermodynamic state";
  output Real der_h "Time derivative of specific enthalpy";

algorithm
  der_h := der_state.h;

  annotation(Inline=true);
end specificEnthalpy_der;
