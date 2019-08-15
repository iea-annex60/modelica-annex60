within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_derp_T
  "Calculates derivative (dh/dp)_T=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dhpT "Specific enthalpy derivative (dh/dp)_T=const.";

protected
  SaturationProperties sat = setSat_T(state.T);
  Real phase_dT = 0;

algorithm
  if state.phase == 0 then
    phase_dT := if not ((state.d < bubbleDensity(sat) and state.d >
                dewDensity(sat)) and state.T <
                fluidConstants[1].criticalTemperature) then 1 else 2;
  else
    phase_dT := state.phase;
  end if;

  if state.phase==1 or phase_dT==1 then
    dhpT := specificEnthalpy_derd_T(state)/pressure_derd_T(state);
  elseif state.phase==2 or phase_dT==2 then
    dhpT := Modelica.Constants.inf;
  end if;

  annotation(Inline=false,
             LateInline=true);
end specificEnthalpy_derp_T;
