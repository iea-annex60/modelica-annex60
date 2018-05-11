within IBPSA.Media.Antifreeze.Validation.BaseClasses;
partial model FluidProperties
  "Partial model that tests the implementation of temperature- and concentration-dependent fluid properties"

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium package";

  parameter Integer n
    "Number of mass fractions to evaluate fluid properties";
  parameter Modelica.SIunits.MassFraction w[n]
    "Mass fraction of additive";
  parameter Modelica.SIunits.Temperature T_min
    "Minimum temperature of mixture";
  parameter Modelica.SIunits.Temperature T_max
    "Maximum temperature of mixture";
  parameter Modelica.SIunits.Temperature reference_T = 293.15
    "Reference temperature";
  Modelica.SIunits.Temperature Tf[n] "Fluid temperature";
  Modelica.SIunits.Density d[n] "Density of fluid mixture";
  Modelica.SIunits.SpecificHeatCapacity cp[n] "Density of fluid mixture";
  Modelica.SIunits.ThermalConductivity lambda[n] "Density of fluid mixture";
  Modelica.SIunits.DynamicViscosity eta[n] "Density of fluid mixture";
  Modelica.SIunits.Temperature T "Temperature";

protected
  parameter Modelica.SIunits.Time dt = 1
    "Simulation length";
  parameter Real convT(unit="K/s") = (T_max-T_min)/dt
    "Rate of temperature change";

equation
  T = T_min + convT*time;
  for i in 1:n loop
    Tf[i] = Medium.polynomialFusionTemperature(w[i],T);
    d[i] = if T >= Tf[i] then Medium.polynomialDensity(w[i],T) else 0.;
    cp[i] = if T >= Tf[i] then Medium.polynomialSpecificHeatCapacityCp(w[i],T) else 0.;
    lambda[i] = if T >= Tf[i] then Medium.polynomialThermalConductivity(w[i],T) else 0.;
    eta[i] = if T >= Tf[i] then Medium.polynomialDynamicViscosity(w[i],T) else 0.;
  end for;

   annotation (
Documentation(info="<html>
<p>
This example checks the implementation of functions that evaluate the
temperature- and concentration-dependent thermophysical properties of the
medium.
</p>
<p>
Thermophysical properties (density, specific heat capacity, thermal conductivity
and dynamic viscosity) are shown as 0 if the temperature is below the fusion
temperature.
</p>
</html>",
revisions="<html>
<ul>
<li>
March 14, 2018, by Massimo Cimmino:<br/>
First implementation.
</li>
</ul>
</html>"));
end FluidProperties;
