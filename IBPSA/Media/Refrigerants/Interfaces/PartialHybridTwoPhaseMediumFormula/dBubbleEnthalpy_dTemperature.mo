within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dBubbleEnthalpy_dTemperature
  "Calculates bubble point enthalpy derivative"
  input SaturationProperties sat "Saturation properties";
  output Real dhlT "Bubble point enthalpy derivative wrt. temperature";

protected
  ThermodynamicState state = setBubbleState(sat);
  Real dpT = 0;

algorithm
  dpT := saturationPressure_derT(sat.Tsat);
  dhlT := specificEnthalpy_derT_p(state) +
          specificEnthalpy_derp_T(state)*dpT;

  annotation(Inline=false,
             LateInline=true);
end dBubbleEnthalpy_dTemperature;
