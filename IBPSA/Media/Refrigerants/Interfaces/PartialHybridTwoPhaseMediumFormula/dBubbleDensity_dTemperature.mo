within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dBubbleDensity_dTemperature
  "Calculates bubble point density derivative"
  input SaturationProperties sat "Saturation properties";
  output Real ddlT "Bubble point density derivative wrt. temperature";

protected
  ThermodynamicState state = setBubbleState(sat);
  Real ddTp = 0;
  Real ddpT = 0;
  Real dpT = 0;

algorithm
  ddTp := density_derT_p(state);
  ddpT := density_derp_T(state);
  dpT := saturationPressure_derT(sat.Tsat);

  ddlT := ddTp + ddpT*dpT;

  annotation(Inline=false,
             LateInline=true);
end dBubbleDensity_dTemperature;
