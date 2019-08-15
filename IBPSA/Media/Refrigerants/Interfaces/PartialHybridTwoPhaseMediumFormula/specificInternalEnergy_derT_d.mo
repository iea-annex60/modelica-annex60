within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificInternalEnergy_derT_d
  "Calculated derivative (du/dT)_d=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real duTd "Derivative (du/dT)_d=const.";

protected
  SaturationProperties sat = setSat_p(state.p);
  Real phase_dT = 0;

  Real delta = 0;
  Real tau = 0;

  Real quality = 0;
  Real dulT = 0;
  Real duvT = 0;
  Real dvlT = 0;
  Real dvvT = 0;
  Real dqualityTv = 0;

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

    duTd := -Modelica.Constants.R/fluidConstants[1].molarMass*
            (tt_fIdg_tt(tau=tau) + tt_fRes_tt(delta=delta,tau=tau));
  else
     quality := (bubbleDensity(sat)/state.d - 1)/
                (bubbleDensity(sat)/dewDensity(sat) - 1);

     dulT := dBubbleInternalEnergy_dTemperature(sat);
     duvT := dDewInternalEnergy_dTemperature(sat);
     dvlT := -1/(bubbleDensity(sat)^2)*dBubbleDensity_dTemperature(sat);
     dvvT := -1/(dewDensity(sat)^2)*dDewDensity_dTemperature(sat);
     dqualityTv := (quality*dvvT + (1-quality)*dvlT)/
                   (1/bubbleDensity(sat)-1/dewDensity(sat));

     duTd := dulT + dqualityTv*((dewEnthalpy(sat)-sat.psat/dewDensity(sat))-
             (bubbleEnthalpy(sat)-sat.psat/bubbleDensity(sat))) +
             quality*(duvT-dulT);
  end if;

  annotation(Inline=false,
             LateInline=true);
end specificInternalEnergy_derT_d;
