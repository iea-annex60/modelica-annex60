within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function density_pT_der
  "Calculates time derivative of density_pT"
  input AbsolutePressure p "Pressure";
  input Temperature T "Specific Enthalpy";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_T "Time derivative of temperature";
  output Real der_d "Time derivative of density";

protected
  ThermodynamicState state = setState_pT(p=p,T=T,phase=phase);

algorithm
  der_d := der_p*density_derp_T(state=state) +
           der_T*density_derT_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end density_pT_der;
