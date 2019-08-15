within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_derp_s
  "Calculates temperature derivative (dT/dp)_s=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dTps "Temperature derivative (dT/dp)_s=const.";

protected
  Real dspT = 0;
  Real dsTp = 0;

algorithm
  dspT := specificEntropy_derd_T(state=state)/pressure_derd_T(state=state);
  dsTp := 1/temperature_ders_p(state=state);
  dTps := -dspT/dsTp;

  annotation(Inline=false,
             LateInline=true);
end temperature_derp_s;
