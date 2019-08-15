within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function jouleThomsonCoefficient
  "Joule-Thomson coefficient of refrigerant"
  input ThermodynamicState state "Thermodynamic state";
  output Real my(unit="K/Pa") "Isothermal throttling coefficient";

algorithm
  my := temperature_derp_h(state);

  annotation(Inline=false,
             LateInline=true);
end jouleThomsonCoefficient;
