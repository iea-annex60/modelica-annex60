within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function setState_psX_der
  "Calculates time derivative of the thermodynamic state record
  calculated by p and s"
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input MassFraction X[:]=reference_X "Mass fractions";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_s "Time derivative of specific entropy";
  input Real der_X[:] "Time derivative of mass fractions";
  output ThermodynamicState der_state
    "Time derivative of thermodynamic state";

algorithm
  der_state := ThermodynamicState(
    d = density_ps_der(p=p,s=s,der_p=der_p,der_s=der_s,phase=phase),
    p = der_p,
    T = temperature_ps_der(p=p,s=s,der_p=der_p,der_s=der_s,phase=phase),
    h = specificEnthalpy_ps_der(p=p,s=s,der_p=der_p,der_s=der_s,phase=phase),
    phase = 0);

  annotation(Inline=true);
end setState_psX_der;
