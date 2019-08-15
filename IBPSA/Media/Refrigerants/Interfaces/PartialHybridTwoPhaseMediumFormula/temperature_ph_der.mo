within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_ph_der
  "Calculates time derivative of temperature_ph"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific enthalpy";
  input FixedPhase phase = 0
    "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_p "Time derivative of pressure";
  input Real der_h "Time derivative of specific enthalpy";
  output Real der_T "Time derivative of density";

protected
  ThermodynamicState state = setState_phX(p=p,h=h,phase=phase);

algorithm
  der_T := der_p*temperature_derp_h(state=state) +
           der_h*temperature_derh_p(state=state);

  annotation(Inline=false,
             LateInline=true);
end temperature_ph_der;
