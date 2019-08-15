within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEnthalpy_derT_d
  "Calculates enthalpy derivative (dh/dT)_d=const"
  input ThermodynamicState state "Thermodynamic state";
  output Real dhTd "Enthalpy derivative (dh/dT)_d=const";

protected
  SaturationProperties sat = setSat_p(state.p);
  Real phase_dT = 0;

  Real delta = 0;
  Real tau = 0;

  Real quality = 0;
  Real dhlT = 0;
  Real dhvT = 0;
  Real dvlT = 0;
  Real dvvT = 0;
  Real dqualityph = 0;

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

    dhTd := Modelica.Constants.R/fluidConstants[1].molarMass*
            (1 - tt_fIdg_tt(tau=tau) - tt_fRes_tt(delta=delta, tau=tau) +
            d_fRes_d(delta=delta,tau=tau) - td_fRes_td(delta=delta, tau=tau));
  elseif state.phase==2 or phase_dT==2 then
    quality := (bubbleDensity(sat)/state.d - 1)/
               (bubbleDensity(sat)/dewDensity(sat) - 1);

    dhlT := dBubbleEnthalpy_dTemperature(sat);
    dhvT := dDewEnthalpy_dTemperature(sat);
    dvlT := -1/(bubbleDensity(sat)^2)*dBubbleDensity_dTemperature(sat);
    dvvT := -1/(dewDensity(sat)^2)*dDewDensity_dTemperature(sat);
    dqualityph := (quality*dvvT + (1-quality)*dvlT)/
                  (1/bubbleDensity(sat)-1/dewDensity(sat));

    dhTd := dhlT + dqualityph*(dewEnthalpy(sat)-bubbleEnthalpy(sat)) +
            quality*(dhvT-dhlT);
  end if;

  annotation(Inline=false,
             LateInline=true);
end specificEnthalpy_derT_d;
