within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function density_derp_s
  "Calculates density derivative (dd/dp)_s=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real ddps "Derivative (dd/dp)_s=const.";

algorithm
    ddps := 1 / (pressure_derd_T(state=state) - pressure_derT_d(state=state)*
            specificEntropy_derd_T(state=state)/
            specificEntropy_derT_d(state=state));

  annotation(Inline=false,
             LateInline=true);
end density_derp_s;
