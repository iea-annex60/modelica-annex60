within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_pT_der
  "Calculates time derivative of specificEnthalpy_pT"
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_T "Time derivative of temperature";
  output Real der_h "Time derivative of specific enthalpy";

protected
  ThermodynamicState state = setState_pT(p=p,T=T,phase=phase);

algorithm
  der_h := der_p*specificEnthalpy_derp_T(state=state) +
           der_T*specificEnthalpy_derT_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end specificEnthalpy_pT_der;
