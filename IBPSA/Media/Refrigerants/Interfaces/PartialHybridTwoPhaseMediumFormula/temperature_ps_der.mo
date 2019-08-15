within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_ps_der
  "Calculates time derivative of temperature_ps"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_s "Time derivative of specific entropy";
  output Real der_T "Time derivative of density";

protected
  ThermodynamicState state = setState_psX(p=p,s=s,phase=phase);

algorithm
  der_T := der_p*temperature_derp_s(state=state) +
           der_s*temperature_ders_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end temperature_ps_der;
