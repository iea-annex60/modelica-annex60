within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function setState_dTX_der
  "Calculates time derivative of the thermodynamic state record
  calculated by d and T"
  input Density d "Density";
  input Temperature T "Temperature";
  input MassFraction X[:]=reference_X "Mass fractions";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_d "Time derivative of density";
  input Real der_T "Time derivative of temperature";
  input Real der_X[:] "Time derivative of mass fractions";
  output ThermodynamicState der_state
    "Time derivative of thermodynamic state";

algorithm
  der_state := ThermodynamicState(
    d = der_d,
    p = pressure_dT_der(d=d,T=T,der_d=der_d,der_T=der_T,phase=phase),
    T = der_T,
    h = specificEnthalpy_dT_der(d=d,T=T,der_d=der_d,der_T=der_T,phase=phase),
    phase = 0);

  annotation(Inline=true);
end setState_dTX_der;
