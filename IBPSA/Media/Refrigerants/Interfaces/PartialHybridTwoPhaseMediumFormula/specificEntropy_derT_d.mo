within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEntropy_derT_d
  "Calculates derivative (ds/dT)_d=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dsTd "Derivative (ds/dT)_d=const.";

protected
    Real delta = state.d/(fluidConstants[1].molarMass/
                 fluidConstants[1].criticalMolarVolume);
    Real tau = fluidConstants[1].criticalTemperature/state.T;

algorithm
  dsTd := 1/state.T*specificHeatCapacityCv(state);

  annotation(Inline=false,
             LateInline=true);
end specificEntropy_derT_d;
