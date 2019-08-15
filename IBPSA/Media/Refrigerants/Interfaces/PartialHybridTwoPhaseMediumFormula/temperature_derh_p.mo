within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_derh_p
  "Calculates temperature derivative (dT/dh)_p=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dThp "Temperature derivative (dT/dh)_p=const.";

algorithm
  dThp := 1/specificEnthalpy_derT_p(state);

  annotation(Inline=false,
             LateInline=true);
end temperature_derh_p;
