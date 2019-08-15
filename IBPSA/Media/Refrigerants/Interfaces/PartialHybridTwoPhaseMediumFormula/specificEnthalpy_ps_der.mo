within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_ps_der
  "Calculates time derivative of specificEnthalpy_ps"
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_s "Time derivative of specific entropy";
  output Real der_h "Time derivative of specific enthalpy";

protected
  ThermodynamicState state = setState_psX(p=p,s=s,phase=phase);

algorithm
  der_h := der_p*specificEnthalpy_derp_s(state=state) +
           der_s*specificEnthalpy_ders_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end specificEnthalpy_ps_der;
