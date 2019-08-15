within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function pressure_derT_d
  "Calculates pressure derivative (dp/dT)_d=const"
  input ThermodynamicState state "Thermodynamic state";
  output Real dpTd "Pressure derivative (dp/dT)_d=const";

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

    dpTd := Modelica.Constants.R/fluidConstants[1].molarMass*state.d*
            (1+d_fRes_d(delta=delta,tau=tau)-td_fRes_td(delta=delta,tau=tau));
  elseif state.phase==2 or phase_dT==2 then
    dpTd := saturationPressure_derT(sat.Tsat);
  end if;

  annotation(Inline=false,
             LateInline=true);
end pressure_derT_d;
