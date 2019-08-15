within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_ders_p
  "Calculates derivative (dh/ds)_p=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dhsp "Specific enthalpy derivative (dh/dT)_p=const.";

algorithm
  dhsp := state.T;

  annotation(Inline=true);
end specificEnthalpy_ders_p;
