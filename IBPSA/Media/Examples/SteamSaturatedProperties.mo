within IBPSA.Media.Examples;
model SteamSaturatedProperties
  "Model that tests the implementation of the steam properties at saturated liquid and vapor states"
  extends Modelica.Icons.Example;

  package Medium = IBPSA.Media.Steam  "Steam medium model";

  parameter Modelica.SIunits.Temperature TMin = 273.16
    "Minimum temperature for the simulation";
  parameter Modelica.SIunits.Temperature TMax = 643.15
    "Maximum temperature for the simulation";
  parameter Modelica.SIunits.Pressure pMin = 620
    "Minimum pressure for the simulation";
  parameter Modelica.SIunits.Pressure pMax = 20000000
    "Maximum pressure for the simulation";
  Modelica.SIunits.Pressure pSat "Saturation pressure";
  Medium.Temperature TSat  "Saturation temperature";
  Modelica.SIunits.Conversions.NonSIunits.Temperature_degC TSat_degC
    "Celsius saturation temperature";

  Medium.SaturationProperties sat "Medium saturation state";

  Modelica.SIunits.Density dl "Density of saturated liquid";
  Modelica.SIunits.Density dv "Density of saturated vapor";
  Modelica.SIunits.SpecificEnthalpy hl "Enthalpy of saturated liquid";
  Modelica.SIunits.SpecificEnthalpy hv "Enthalpy of saturated vapor";
  Modelica.SIunits.SpecificEnthalpy hlv "Enthalpy of vaporization";
  Modelica.SIunits.SpecificEnthalpy hlv_old "Enthalpy of vaporization";
  Modelica.SIunits.SpecificEntropy sl "Entropy of saturated liquid";
  Modelica.SIunits.SpecificEntropy sv "Entropy of saturated vapor";
  Modelica.SIunits.SpecificEntropy slv "Entropy of vaporization";

protected
  constant Real conv(unit="1/s") = 1 "Conversion factor to satisfy unit check";

equation
    // Compute temperatures that are used as input to the functions
    pSat = pMin + conv*time * (pMax-pMin);
    TSat_degC = Modelica.SIunits.Conversions.to_degC(TSat);

    // Saturation state
    sat = Medium.setSat_p(pSat);
    TSat = sat.Tsat;
    assert(TSat > TMin, "Temperature exceeded minimum value.\n" +
      "   TSat = " + String(TSat));
    assert(TSat < TMax, "Temperature exceeded maximum value.\n" +
      "   TSat = " + String(TSat));
    // Check the implementation of the functions
    dl = Medium.densityOfSaturatedLiquid(sat);
    dv = Medium.densityOfSaturatedVapor(sat);
    hl = Medium.enthalpyOfSaturatedLiquid(sat);
    hv = Medium.enthalpyOfSaturatedVapor(sat);
    hlv = Medium.enthalpyOfVaporization(sat);
    hlv_old = Medium.enthalpyOfVaporization_old(sat.Tsat);
    sl = Medium.entropyOfSaturatedLiquid(sat);
    sv = Medium.entropyOfSaturatedVapor(sat);
    slv = Medium.entropyOfVaporization(sat);

   annotation(experiment(Tolerance=1e-6, StopTime=1.0),
__Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Media/Examples/SteamSaturatedProperties.mos"
        "Simulate and plot"),
      Documentation(info="<html>
<p>
This example checks the saturation properties of the medium. 
</p>
<p>
The steam medium is designed for single phase (saturated or 
superheated) vapor (x=1). The saturated state functions are 
provided for models involving phase change by implementing both
the <a href=\"modelica://IBPSA.Media.Steam\">
IBPSA.Media.Steam</a> model (vapor phase) and the
<a href=\"modelica://IBPSA.Media.Water\"> IBPSA.Media.Water</a> 
model (liquid phase). See <a href=\"modelica://IBPSA.Media.Steam\">
IBPSA.Media.Steam</a> for more information.
</p>
</html>",
revisions="<html>
<ul>
<li>
March 4, 2020, by Kathryn Hinkelman:<br/>
First implementation.
</li>
</ul>
</html>"));
end SteamSaturatedProperties;