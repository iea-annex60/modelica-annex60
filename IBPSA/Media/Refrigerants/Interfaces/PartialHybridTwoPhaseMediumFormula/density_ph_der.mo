within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function density_ph_der
  "Calculates time derivative of density_ph"
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific Enthalpy";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_h "Time derivative of specific enthalpy";
  output Real der_d "Time derivative of density";

protected
  ThermodynamicState state = setState_phX(p=p,h=h,phase=phase);

algorithm
  der_d := der_p*density_derp_h(state=state) +
           der_h*density_derh_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end density_ph_der;
