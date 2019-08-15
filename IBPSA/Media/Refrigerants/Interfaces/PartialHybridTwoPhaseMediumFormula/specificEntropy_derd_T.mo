within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function specificEntropy_derd_T
  "Calculates derivative (ds/dd)_T=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real dsdT "Derivative (ds/dd)_T=const.";

protected
  SaturationProperties sat = setSat_p(state.p);
  Real phase_dT = 0;

  Real R = Modelica.Constants.R/fluidConstants[1].molarMass;
  Real delta = 0;
  Real tau = 0;

  ThermodynamicState state_l;
  ThermodynamicState state_v;
  Real delta_l = 0;
  Real tau_l = 0;
  Real delta_v = 0;
  Real tau_v = 0;
  Real dslTd = 0;
  Real dsvTd = 0;
  Real dslpT = 0;
  Real dsvpT = 0;
  Real dslTp = 0;
  Real dsvTp = 0;
  Real dsldT = 0;
  Real dsvdT = 0;
  Real dpT = 0;
  Real dsvT = 0;
  Real dslT = 0;
  Real quality = 0;
  Real dvT_liq = 0;
  Real dvT_vap = 0;
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

    dsdT := -R/state.d*(1+d_fRes_d(tau=tau, delta=delta) -
            td_fRes_td(tau=tau, delta=delta));
  else
     state_l := setBubbleState(sat);
     state_v := setDewState(sat);
     delta_l := state_l.d/(fluidConstants[1].molarMass/
                fluidConstants[1].criticalMolarVolume);
     delta_v := state_v.d/(fluidConstants[1].molarMass/
                fluidConstants[1].criticalMolarVolume);
     tau_l := fluidConstants[1].criticalTemperature/state_l.T;
     tau_v := fluidConstants[1].criticalTemperature/state_v.T;

     dslTd := specificEntropy_derd_T(state_l);
     dsvTd := specificEntropy_derd_T(state_v);
     dsldT := R/state_l.d*(-(1 + d_fRes_d(tau=tau_l,delta=delta_l)) +
              (td_fRes_td(tau=tau_l,delta=delta_l)));
     dsvdT := R/state_v.d*(-(1 + d_fRes_d(tau=tau_v,delta=delta_v)) +
              (td_fRes_td(tau=tau_v,delta=delta_v)));
     dslTp := dslTd - dsldT*pressure_derT_d(state_l)/pressure_derd_T(state_l);
     dsvTp := dsvTd - dsvdT*pressure_derT_d(state_v)/pressure_derd_T(state_v);
     dslpT := dslpT/pressure_derd_T(state_l);
     dsvpT := dsvpT/pressure_derd_T(state_v);
     dpT := saturationPressure_derT(sat.Tsat);

     dsvT := dsvTp + dsvpT*dpT;
     dslT := dslTp + dslpT*dpT;

     quality := (bubbleDensity(sat)/state.d-1)/
                (bubbleDensity(sat)/dewDensity(sat)-1);
     dvT_liq := -1/state_l.d^2*(density_derT_p(state_l) +
                density_derp_T(state_l)*dpT);
     dvT_vap := -1/state_v.d^2*(density_derT_p(state_v) +
                density_derp_T(state_v)*dpT);
     dqualityTv := (-dvT_liq*(1/state_v.d - 1/state_l.d) -
                   (1/state.d - 1/state_l.d)*(dvT_vap - dvT_liq))/
                   (1/state_v.d - 1/state_l.d)^2;

     dsdT := dslT + quality*(dsvT - dslT) + dqualityTv*
             (dewEntropy(sat)-bubbleEntropy(sat));
  end if;

  annotation(Inline=false,
             LateInline=true);
end specificEntropy_derd_T;
