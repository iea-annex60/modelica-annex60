within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function isothermalThrottlingCoefficient
  "Isothermal throttling coefficient of refrigerant"
  input ThermodynamicState state "Thermodynamic state";
  output Real delta_T(unit="J/(Pa.kg)") "Isothermal throttling coefficient";

algorithm
  delta_T := specificEnthalpy_derp_T(state);

  annotation(Inline=false,
             LateInline=true);
end isothermalThrottlingCoefficient;
