within Annex60.Airflow.Multizone;
model MediumColumn
  "Vertical shaft with no friction and no storage of heat and mass"
  import Modelica.Constants;

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium in the component" annotation (choicesAllMatching=true);

  parameter SI.Length h(min=0) = 3 "Height of shaft";
  parameter Annex60.Airflow.Multizone.Types.densitySelection densitySelection
    "Select how to pick density" annotation (Evaluate=true);
  parameter Boolean allowFlowReversal=true
    "= true to allow flow reversal, false restricts to design direction (port_a -> port_b)"
    annotation (Dialog(tab="Assumptions"),Evaluate=true);

  Modelica.Fluid.Interfaces.FluidPort_a port_a(
    redeclare package Medium = Medium,
    m_flow(min=if allowFlowReversal then -Constants.inf else 0),
    p(start=Medium.p_default))
    "Fluid connector a (positive design flow direction is from port_a to port_b)"
    annotation (Placement(transformation(extent={{-10,90},{10,110}}),
        iconTransformation(extent={{-10,90},{10,110}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_b(
    redeclare package Medium = Medium,
    m_flow(max=if allowFlowReversal then +Constants.inf else 0),
    p(start=Medium.p_default))
    "Fluid connector b (positive design flow direction is from port_a to port_b)"
    annotation (Placement(transformation(extent={{10,-110},{-10,-90}}), iconTransformation(extent={{10,-110},{-10,-90}})));

  SI.VolumeFlowRate V_flow
    "Volume flow rate at inflowing port (positive when flow from port_a to port_b)";
  SI.MassFlowRate m_flow(start=0)
    "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
  SI.Pressure dp(start=0, displayUnit="Pa")
    "Pressure difference between port_a and port_b";
  SI.Density rho "Density in medium column";
protected
  Medium.ThermodynamicState sta_a=Medium.setState_phX(
      port_a.p,
      actualStream(port_a.h_outflow),
      actualStream(port_a.Xi_outflow)) "Medium properties in port_a";
  Medium.MassFraction Xi[Medium.nXi] "Mass fraction used to compute density";
initial equation
  /*
   assert(abs(Medium.density(Medium.setState_pTX(
    Medium.p_default,
    Medium.T_default,
    Medium.X_default)) - Medium.density(Medium.setState_pTX(
    Medium.p_default,
    Medium.T_default + 5,
    Medium.X_default))) > 1E-10,
    "Error: The density of the medium that is used to compute buoyancy force is independent of temperature."
    + "\n       You need to select a different medium model.");
    */
  // The next assert tests for all allowed values of the enumeration.
  // Testing against densitySelection > 0 gives an error in OpenModelica as enumerations start with 1.
  assert(densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.fromTop
     or densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.fromBottom
     or densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.actual,
    "You need to set the parameter \"densitySelection\" for the model MediumColumn.");
equation
  // Design direction of mass flow rate
  m_flow = port_a.m_flow;

  // Pressure difference between ports
  // Xi is computed first as it is used in two expression, and in one
  // of them only one component is used.
  // We test for Medium.nXi == 0 as Modelica.Media.Air.SimpleAir has no
  // moisture and hence Xi[1] is an illegal statement.
  // We first compute temperature and then invoke a density function that
  // takes temperature as an argument. Simply calling a density function
  // of a medium that takes enthalpy as an argument would be dangerous
  // as different media can have different datum for the enthalpy.
  if (densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.fromTop) then
      Xi = inStream(port_a.Xi_outflow);
      rho = Annex60.Utilities.Psychrometrics.Functions.density_pTX(
        p=port_a.p,
        T=Medium.temperature(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), Xi)),
        X_w=if Medium.nXi == 0 then 0 else Xi[1]);
  elseif (densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.fromBottom) then
      Xi = inStream(port_b.Xi_outflow);
      rho = Annex60.Utilities.Psychrometrics.Functions.density_pTX(
        p=port_b.p,
        T=Medium.temperature(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), Xi)),
        X_w=if Medium.nXi == 0 then 0 else Xi[1]);
   else
      Xi = actualStream(port_a.Xi_outflow);
      rho = Annex60.Utilities.Psychrometrics.Functions.density_pTX(
        p=port_a.p,
        T=Medium.temperature(Medium.setState_phX(port_a.p, actualStream(port_a.h_outflow), Xi)),
        X_w=if Medium.nXi == 0 then 0 else Xi[1]);
  end if;

  V_flow = m_flow/Medium.density(sta_a);

  dp = port_a.p - port_b.p;
  dp = -h*rho*Modelica.Constants.g_n;

  // Isenthalpic state transformation (no storage and no loss of energy)
  port_a.h_outflow = inStream(port_b.h_outflow);
  port_b.h_outflow = inStream(port_a.h_outflow);

  // Mass balance (no storage)
  port_a.m_flow + port_b.m_flow = 0;

  // Transport of substances
  port_a.Xi_outflow = inStream(port_b.Xi_outflow);
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);

  port_a.C_outflow = inStream(port_b.C_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);

  annotation (
    Icon(graphics={
        Line(
          points={{0,100},{0,-100},{0,-98}}),
        Text(
          extent={{24,-78},{106,-100}},
          lineColor={0,0,127},
          textString="Bottom"),
        Text(
          extent={{32,104},{98,70}},
          lineColor={0,0,127},
          textString="Top"),
        Text(
          extent={{36,26},{88,-10}},
          lineColor={0,0,127},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid,
          textString="h=%h"),
        Rectangle(
          extent={{-16,80},{16,-80}},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Rectangle(
          visible=densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.fromTop,
          extent={{-16,80},{16,0}},
          fillColor={85,170,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Text(
          extent={{-50.5,20.5},{50.5,-20.5}},
          lineColor={0,0,127},
          origin={-72.5,-12.5},
          rotation=90,
          textString="%name"),
        Rectangle(
          visible=densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.actual,
          extent={{-16,80},{16,54}},
          fillColor={85,170,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Rectangle(
          visible=densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.fromBottom,
          extent={{-16,0},{16,-82}},
          fillColor={85,170,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Rectangle(
          visible=densitySelection == Annex60.Airflow.Multizone.Types.densitySelection.actual,
          extent={{-16,-55},{16,-80}},
          fillColor={85,170,255},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None,
          lineColor={0,0,0})}),
defaultComponentName="col",
Documentation(info="<html>
<p>
This model describes the pressure difference of a vertical medium
column. It can be used to model the pressure difference caused by
stack effect.
</p>
<p>
The model can be used with the following three configurations, which are
controlled by the setting of the parameter <code>densitySelection</code>:
</p>
<ul>
<li>
<code>top</code>:
Use this setting to use the density from the volume that is connected
to <code>port_a</code>.
</li>
<li>
<code>bottom</code>:
Use this setting to use the density from the volume that is connected
to <code>port_b</code>.
</li>
<li>
<code>actual</code>:
Use this setting to use the density based on the actual flow direction.
</li>
</ul>
<p>
The settings <code>top</code> and <code>bottom</code>
should be used when rooms or different floors of a building are
connected since multizone airflow models assume that each floor is completely mixed.
For these two seetings, this model will compute the pressure between the center of the room
and an opening that is at height <code>h</code> relative to the center of the room.
The setting <code>actual</code> may be used to model a chimney in which
a column of air will change its density based on the flow direction.
</p>
<p>
In this model, the parameter <code>h</code> must always be positive, and the port <code>port_a</code> must be
at the top of the column.
</p>
<p>
For a steady-state model, use
<a href=\"modelica://Annex60.Airflow.Multizone.MediumColumnDynamic\">
Annex60.Airflow.Multizone.MediumColumnDynamic</a> instead of this model.
</p>
</html>",
revisions="<html>
<ul>
<li>
February 24, 2015 by Michael Wetter:<br/>
Changed model to use
<a href=\"modelica://Annex60.Utilities.Psychrometrics.Functions.density_pTX\">
Annex60.Utilities.Psychrometrics.Functions.density_pTX</a>
for the density computation
as
<a href=\"modelica://Annex60.Media.Air.density\">
Annex60.Media.Air.density</a>
does not depend on temperature.
</li>
<li>
December 22, 2014 by Michael Wetter:<br/>
Removed <code>Modelica.Fluid.System</code>
to address issue
<a href=\"https://github.com/lbl-srg/modelica-Annex60/issues/311\">#311</a>.
</li>
<li><i>October 4, 2014</i> by Michael Wetter:<br/>
Removed assignment of <code>port_?.p.nominal</code> to avoid a warning in OpenModelica because
alias sets have different nominal values.
</li>
<li><i>April 17, 2013</i> by Michael Wetter:<br/>
       Reformulated the assert statement that checks for the correct value of <code>densitySelection</code>.
</li>
<li><i>July 28, 2010</i> by Michael Wetter:<br/>
       Changed sign for pressure difference.
</li>
<li><i>July 20, 2010</i> by Michael Wetter:<br/>
       Migrated model to Modelica 3.1 and integrated it into the Buildings library.
       Reimplemented assignment of density based on flow direction or based on outflowing state.
</li>
<li><i>February 24, 2005</i> by Michael Wetter:<br/>
       Released first version.
</ul>
</html>"));
end MediumColumn;
