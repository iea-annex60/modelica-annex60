within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_ders_p
  "Calculates temperature derivative (dT/ds)@p=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dTsp "Temperature derivative (dT/ds)@p=const.";

algorithm
  dTsp := 1/(specificEntropy_derT_d(state)-specificEntropy_derd_T(state)*
          pressure_derT_d(state)/pressure_derd_T(state));

  annotation(Inline=false,
             LateInline=true);
end temperature_ders_p;
