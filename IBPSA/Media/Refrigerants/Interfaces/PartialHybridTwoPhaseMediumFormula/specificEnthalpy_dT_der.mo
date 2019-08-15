within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_dT_der
  "Calculates time derivative of specificEnthalpy_dT"
  input Density d "Density";
  input Temperature T "Temperature";
  input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_d "Time derivative of density";
  input Real der_T "Time derivative of temperature";
  output Real der_h "Time derivative of specific enthalpy";

protected
  ThermodynamicState state = setState_dT(d=d,T=T,phase=phase);

algorithm
  der_h := der_d*specificEnthalpy_derd_T(state=state) +
           der_T*specificEnthalpy_derT_d(state=state);

  annotation(Inline=false,
             LateInline=true);
end specificEnthalpy_dT_der;
