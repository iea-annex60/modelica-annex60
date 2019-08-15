within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function setState_pTX_der
  "Calculates time derivative of the thermodynamic state record
  calculated by p and T"
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input MassFraction X[:]=reference_X "Mass fractions";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_T "Time derivative of temperature";
  input Real der_X[:] "Time derivative of mass fractions";
  output ThermodynamicState der_state
    "Time derivative of thermodynamic state";

algorithm
  der_state := ThermodynamicState(
    d = density_pT_der(p=p,T=T,der_p=der_p,der_T=der_T,phase=phase),
    p = der_p,
    T = der_T,
    h = specificEnthalpy_pT_der(p=p,T=T,der_p=der_p,der_T=der_T,phase=phase),
    phase = 0);

  annotation(Inline=true);
end setState_pTX_der;
