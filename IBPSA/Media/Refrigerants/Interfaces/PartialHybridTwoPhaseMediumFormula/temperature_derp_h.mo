within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function temperature_derp_h
  "Calculates temperature derivative (dT/dp)_h=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dTph "Temperature derivative (dT/dp)_h=const.";

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
    dTph := 1 / (pressure_derT_d(state) - pressure_derd_T(state) *
            specificEnthalpy_derT_d(state)/specificEnthalpy_derd_T(state));
  elseif state.phase==2 or phase_dT==2 then
    dTph := saturationTemperature_derp(sat.psat);
  end if;

  annotation(Inline=false,
             LateInline=true);
end temperature_derp_h;
