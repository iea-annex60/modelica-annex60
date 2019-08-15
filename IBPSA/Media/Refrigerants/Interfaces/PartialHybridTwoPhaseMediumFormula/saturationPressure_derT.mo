within IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula;
function saturationPressure_derT
  "Calculates derivative (dp/dT)_saturation"
  input Temperature T "Temperature";
  output Real dpT "Saturation temperature derivative (dp/dT)_saturation";

protected
  SaturationProperties sat = setSat_T(T=T);

algorithm
  dpT := (dewEntropy(sat)-bubbleEntropy(sat))/
         (1/dewDensity(sat)-1/bubbleDensity(sat));

  annotation(Inline=true);
end saturationPressure_derT;
