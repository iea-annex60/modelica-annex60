within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function setState_phX_der
  "Calculates time derivative of the thermodynamic state record
  calculated by p and h"
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific enthalpy";
  input MassFraction X[:]=reference_X "Mass fractions";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_h "Time derivative of specific enthalpy";
  input Real der_X[:] "Time derivative of mass fractions";
  output ThermodynamicState der_state
    "Time derivative of thermodynamic state";

algorithm
  der_state := ThermodynamicState(
    d = density_ph_der(p=p,h=h,der_p=der_p,der_h=der_h,phase=phase),
    p = der_p,
    T = temperature_ph_der(p=p,h=h,der_p=der_p,der_h=der_h,phase=phase),
    h = der_h,
    phase = 0);

  annotation(Inline=true);
end setState_phX_der;
