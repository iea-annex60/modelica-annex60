within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dDewEnthalpy_dTemperature
  "Calculates dew point enthalpy derivative"
  input SaturationProperties sat "Saturation properties";
  output Real dhvT "Dew point enthalpy derivative wrt. temperature";

protected
  ThermodynamicState state = setDewState(sat);
  Real dpT = 0;

algorithm
  dpT := saturationPressure_derT(sat.Tsat);
  dhvT := specificEnthalpy_derT_p(state) +
          specificEnthalpy_derp_T(state)*dpT;

  annotation(Inline=false,
             LateInline=true);
end dDewEnthalpy_dTemperature;
