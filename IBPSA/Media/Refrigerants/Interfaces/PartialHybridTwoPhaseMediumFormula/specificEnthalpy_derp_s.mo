within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_derp_s
  "Calculates derivative (dh/dp)_s=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dhps "Specific enthalpy derivative (dh/dp)_s=const.";

algorithm
  dhps := 1/state.d;

  annotation(Inline=true);
end specificEnthalpy_derp_s;
