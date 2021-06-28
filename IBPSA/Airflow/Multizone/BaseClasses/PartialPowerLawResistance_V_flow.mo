within IBPSA.Airflow.Multizone.BaseClasses;
partial model PartialPowerLawResistance_V_flow
  "Flow resistance that uses the power law for computing volumetric flow rate"

  extends IBPSA.Airflow.Multizone.BaseClasses.PartialOneWayFlowElement(
  V_flow = IBPSA.Airflow.Multizone.BaseClasses.powerLawFixedM(
      k=k,
      dp=dp,
      m=m,
      a=a,
      b=b,
      c=c,
      d=d,
      dp_turbulent=dp_turbulent),
    final m_flow_nominal=rho_default*k*dp_turbulent,
    final m_flow_small=1E-4*abs(m_flow_nominal));
   extends IBPSA.Airflow.Multizone.BaseClasses.PowerLawResistanceParameters;



  parameter Real k "Flow coefficient, k = V_flow/ dp^m";

  annotation (
    Documentation(info="<html>
<p>
This model describes the mass flow rate and pressure difference relation
of an orifice in the form
</p>
<pre>
    V_flow = k * dp^m,
</pre>
<p>
where <code>k</code> is a variable and
<code>m</code> a parameter.
For turbulent flow, set <code>m=1/2</code> and
for laminar flow, set <code>m=1</code>.
</p>
<p>
The model is used as a base for the interzonal air flow models.
</p>
</html>",
revisions="<html>
<ul>
<li>
May 12, 2020, by Michael Wetter:<br/>
Changed assignment of <code>m_flow_small</code> to <code>final</code>.
This quantity are not used in this model and models that extend from it.
Hence there is no need for the user to change the value.
</li>
<li>
April 14, 2020, by Michael Wetter:<br/>
Changed <code>homotopyInitialization</code> to a constant.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1341\">IBPSA, #1341</a>.
</li>
<li>
June 24, 2018, by Michael Wetter:<br/>
Removed parameter <code>A</code> because
<a href=\"modelica://IBPSA.Airflow.Multizone.EffectiveAirLeakageArea\">
IBPSA.Airflow.Multizone.EffectiveAirLeakageArea</a>
uses the effective leakage area <code>L</code> rather than <code>A</code>.<br/>
Removed calculation <code>v=V_flow/A</code> as parameter <code>A</code> has been removed.<br/>
Removed parameter <code>lWet</code> as this is only used to compute
the Reynolds number, and the Reynolds number is not used by this model.
Also removed the variable <code>Re</code> for the Reynolds number.<br/>
This change is non-backward compatible.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/932\">IBPSA, #932</a>.
</li>
<li>
May 1, 2018, by Filip Jorissen:<br/>
Set <code>final allowFlowReversal=true</code>.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/877\">#877</a>.
</li>
<li>
January 22, 2016, by Michael Wetter:<br/>
Corrected type declaration of pressure difference.
This is
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/404\">#404</a>.
</li>
<li>
November 19, 2015, by Michael Wetter:<br/>
Removed assignment of parameter
<code>showDesignFlowDirection</code> in <code>extends</code> statement.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/349\">#349</a>.
</li>
<li>
January 21, 2015 by Michael Wetter:<br/>
Changed type of <code>mExc</code> as <code>Modelica.SIunits.Mass</code>
sets <code>min=0</code>, but <code>mExc</code> can be negative.
</li>
<li>
October 8, 2013 by Michael Wetter:<br/>
Changed the parameter <code>useConstantDensity</code> to
<code>useDefaultProperties</code> and also applied the parameter
to the computation of the dynamic viscosity.
The conversion script can be used to update this parameter.<br/>
Change model to not use the instance <code>sta_a</code>, as this
may be conditionally removed and hence it is not proper Modelica
syntax to use it outside of a <code>connect</code> statement.
</li>
<li>
March 27, 2013 by Michael Wetter:<br/>
Added assignment of initial value for <code>mExc</code> to avoid error when checking model
in pedantic mode with Dymola 2014.
</li>
<li>
December 14, 2012 by Michael Wetter:<br/>
Renamed protected parameters for consistency with the naming conventions.
</li>
<li>December 6, 2011 by Michael Wetter:<br/>
       Removed <code>fixed=false</code> attribute of protected parameter
       <code>k</code>.
</li>
<li>July 20, 2010 by Michael Wetter:<br/>
       Migrated model to Modelica 3.1 and integrated it into the Buildings library.
</li>
<li>February 4, 2005 by Michael Wetter:<br/>
       Released first version.
</ul>
</html>"));
end PartialPowerLawResistance_V_flow;
