within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificInternalEnergy_derd_T
  "Calculated derivative (du/dd)_T=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dudT "Derivative (du/dd)_T=const.";

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

  if state.phase == 1 or phase_dT == 1 then
    delta := state.d/(fluidConstants[1].molarMass/
             fluidConstants[1].criticalMolarVolume);
    tau := fluidConstants[1].criticalTemperature/state.T;
    dudT := Modelica.Constants.R/fluidConstants[1].molarMass*state.T/state.d
            *td_fRes_td(tau=tau, delta=delta);
  else
    dudT := -1/state.d^2*((bubbleEnthalpy(sat)-sat.psat/bubbleDensity(sat))-
            (dewEnthalpy(sat)-sat.psat/dewDensity(sat)))/
            (1/bubbleDensity(sat)-1/dewDensity(sat));
  end if;

  annotation(Inline=false,
             LateInline=true);
end specificInternalEnergy_derd_T;
