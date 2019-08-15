within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function density_ders_p
  "Calculates density derivative (dd/ds)_p=const."
  input ThermodynamicState state "Thermodynamic state";
  output Real ddsp "Derivative (dd/ds)_p=const.";

algorithm
    ddsp := density_derh_p(state=state)*state.T;

  annotation(Inline=false,
             LateInline=true);
end density_ders_p;
