within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dBubbleInternalEnergy_dTemperature
  "Calculates bubble point internal energy derivative"
  input SaturationProperties sat "Saturation properties";
  output Real dulT
    "Bubble point internal energy derivative wrt. temperature";

protected
  ThermodynamicState state = setBubbleState(sat);
  Real duTd = 0;
  Real dudT = 0;
  Real duTp = 0;
  Real dupT = 0;
  Real dpT = 0;

algorithm
  duTd := specificInternalEnergy_derT_d(state);
  dudT := specificInternalEnergy_derd_T(state);
  duTp := duTd - dudT*pressure_derT_d(state)/pressure_derd_T(state);
  dupT := dudT/pressure_derd_T(state);
  dpT := saturationPressure_derT(sat.Tsat);

  dulT := duTp +dupT*dpT;

  annotation(Inline=false,
             LateInline=true);
end dBubbleInternalEnergy_dTemperature;
