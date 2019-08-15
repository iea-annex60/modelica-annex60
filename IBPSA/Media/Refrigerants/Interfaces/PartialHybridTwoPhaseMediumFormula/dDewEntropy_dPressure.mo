within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function dDewEntropy_dPressure
  "Calculates dew point entropy derivative"
  input SaturationProperties sat "Saturation properties";
  output Real dsvdp "Dew point entropy derivative at constant Temperature";

protected
  ThermodynamicState state = setDewState(sat);
  Real dsdT = 0;
  Real dpdT = 0;
  Real dsTd = 0;
  Real dpTd = 0;
  Real dsTp = 0;
  Real dspT = 0;
  Real dTp = 0;

algorithm
  dsdT := specificEntropy_derd_T(state);
  dsTd := specificEntropy_derT_d(state);
  dsTp := dsTd - dsdT*pressure_derT_d(state)/pressure_derd_T(state);
  dspT := dsdT/pressure_derd_T(state);
  dTp := saturationTemperature_derp(sat.psat);

  dsvdp := dspT + dsTp * dTp;

  annotation(Inline=false,
             LateInline=true);
end dDewEntropy_dPressure;
