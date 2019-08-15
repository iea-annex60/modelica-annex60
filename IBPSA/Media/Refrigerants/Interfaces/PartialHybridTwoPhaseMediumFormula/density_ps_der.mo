within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function density_ps_der
  "Calculates time derivative of density_ps"
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific Entropy";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_s "Time derivative of specific entropy";
  output Real der_d "Time derivative of density";

protected
  ThermodynamicState state = setState_ps(p=p,s=s,phase=phase);

algorithm
  der_d := der_p*density_derp_s(state=state) +
           der_s*density_ders_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end density_ps_der;
