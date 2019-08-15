within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_derT_p
  "Calculates derivative (dh/dT)_p=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dhTp "Specific enthalpy derivative (dh/dT)_p=const.";

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
    dhTp := specificEnthalpy_derT_d(state)-specificEnthalpy_derd_T(state)*
            pressure_derT_d(state)/pressure_derd_T(state);
  elseif state.phase==2 or phase_dT==2 then
    dhTp := Modelica.Constants.inf;
  end if;

  annotation(Inline=false,
             LateInline=true);
end specificEnthalpy_derT_p;
