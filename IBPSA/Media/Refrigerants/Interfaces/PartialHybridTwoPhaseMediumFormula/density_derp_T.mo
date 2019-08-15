within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function extends density_derp_T
  "Calculates the derivative (dd/dp)_T=const."
protected
  SaturationProperties sat = setSat_p(state.p);
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
    ddpT := 1/pressure_derd_T(state);
  elseif state.phase==2 or phase_dT==2 then
    ddpT := Modelica.Constants.inf;
  end if;

  annotation(Inline=false,
             LateInline=true);
end density_derp_T;
