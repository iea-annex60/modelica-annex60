within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function pressure_derd_T
  "Calculates pressure derivative (dp/dd)_T=const"
  input ThermodynamicState state "Thermodynamic state";
  output Real dpdT "Pressure derivative (dp/dd)_T=const";

protected
  SaturationProperties sat = setSat_p(state.p);
  Real phase_dT = 0;

  Real delta = 0;
  Real tau = 0;

algorithm
  if state.phase == 0 then
    phase_dT := if not ((state.d < bubbleDensity(sat) and state.d >
                dewDensity(sat)) and state.T <
                fluidConstants[1].criticalTemperature) then 1 else 2;
  else
    phase_dT := state.phase;
  end if;

  if state.phase==1 or phase_dT==1 then
    delta := state.d/(fluidConstants[1].molarMass/
             fluidConstants[1].criticalMolarVolume);
    tau := fluidConstants[1].criticalTemperature/state.T;

    dpdT := Modelica.Constants.R/fluidConstants[1].molarMass*
            state.T*(1 + 2*d_fRes_d(delta=delta,tau=tau) +
            dd_fRes_dd(delta=delta,tau=tau));
  elseif state.phase==2 or phase_dT==2 then
    dpdT := Modelica.Constants.small;
  end if;

  annotation(Inline=true);
end pressure_derd_T;
