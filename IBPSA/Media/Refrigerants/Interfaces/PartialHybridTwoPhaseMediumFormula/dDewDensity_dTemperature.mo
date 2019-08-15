within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dDewDensity_dTemperature
  "Calculates dew point debsity derivative"
  input SaturationProperties sat "Saturation properties";
  output Real ddvT "Dew point density derivative wrt. temperature";

protected
  ThermodynamicState state = setDewState(sat);
  Real ddTp = 0;
  Real ddpT = 0;
  Real dpT = 0;

algorithm
  ddTp := density_derT_p(state);
  ddpT := density_derp_T(state);
  dpT := saturationPressure_derT(sat.Tsat);

  ddvT := ddTp + ddpT*dpT;

  annotation(Inline=false,
             LateInline=true);
end dDewDensity_dTemperature;
