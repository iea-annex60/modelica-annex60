within IBPSA.Airflow.Multizone.Examples;
model PressurizationData
  "Model showing how the 'Powerlaw_1DataPoint' model can be used when data is available from a pressurization test."
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Air;


  parameter Real n50=3 "ACH50, air changes at 50Pa";

  BoundaryConditions.WeatherData.ReaderTMY3       weaDat(filNam=
        Modelica.Utilities.Files.loadResource("modelica://IBPSA/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos"))
    annotation (Placement(transformation(extent={{-96,6},{-76,26}})));
  Fluid.Sources.Outside_CpLowRise       west(
    redeclare package Medium = Medium,
    s=5,
    azi=IBPSA.Types.Azimuth.W,
    Cp0=0.6,
    nPorts=1) "Model with outside conditions"
    annotation (Placement(transformation(extent={{84,6},{64,26}})));
  Fluid.Sources.Outside_CpLowRise east(
    redeclare package Medium = Medium,
    s=5,
    azi=IBPSA.Types.Azimuth.E,
    Cp0=0.6,
    nPorts=1) "Model with outside conditions"
    annotation (Placement(transformation(extent={{-64,6},{-44,26}})));
  Fluid.MixingVolumes.MixingVolume       Room(
    redeclare package Medium = Medium,
    V=2.5*5*5,
    nPorts=2,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    m_flow_nominal=0.01)
    annotation (Placement(transformation(extent={{0,20},{20,40}})));
  Powerlaw_1Datapoint powlaw_1dat(redeclare package Medium = Medium,
    m=0.66,                                                          dP1(
        displayUnit="Pa") = 50,
    m1_flow=0.5*(Room.V*n50*1.2))
    "Crack in envelope representing 50% of the leakage area"
    annotation (Placement(transformation(extent={{-34,6},{-14,26}})));
  Powerlaw_1Datapoint
                   powlaw_1dat1(
                            redeclare package Medium = Medium,
    m=0.66,
    dP1(displayUnit="Pa") = 50,
    m1_flow=0.5*(Room.V*n50*1.2))
    "Crack in envelope representing 50% of the leakage area"
    annotation (Placement(transformation(extent={{34,6},{54,26}})));
equation
  connect(weaDat.weaBus, west.weaBus) annotation (Line(
      points={{-76,16},{-68,16},{-68,-18},{84,-18},{84,16.2}},
      color={255,204,51},
      thickness=0.5));
  connect(east.weaBus, weaDat.weaBus) annotation (Line(
      points={{-64,16.2},{-64,16},{-76,16}},
      color={255,204,51},
      thickness=0.5));
  connect(east.ports[1], powlaw_1dat.port_a)
    annotation (Line(points={{-44,16},{-34,16}},   color={0,127,255}));
  connect(powlaw_1dat.port_b, Room.ports[1])
    annotation (Line(points={{-14,16},{8,16},{8,20}},      color={0,127,255}));
  connect(powlaw_1dat1.port_a, Room.ports[2])
    annotation (Line(points={{34,16},{12,16},{12,20}},    color={0,127,255}));
  connect(powlaw_1dat1.port_b, west.ports[1])
    annotation (Line(points={{54,16},{64,16}},   color={0,127,255}));
  annotation (__Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/Examples/OneEffectiveAirLeakageArea.mos"
        "Simulate and plot"),
        experiment(
      StopTime=2592000,
      Interval=600,
      Tolerance=1e-08,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>
This model illustrates the use of the Powerlaw_1DataPoint to model infiltration through the building evelope for a known <i>n50 value</i> (also known as ACH50).
As the n50 value and the building volume is known, the flow at 50Pa is known. Dividing this flow accross the entire envelope (typically surface weighted) and using this Powerlaw_1DataPoint, the infiltration airflow at lower pressure differences can be modelled.
<br>
In this example, the two models each represent 50% of the surface where airflow occured due to the pressurization test.

</p>
</html>", revisions="<html>
<ul>
<li>
May 03, 2021 by Klaas De Jonge:<br/>
Added example for simulating infiltration airflow using the Powerlaw_1DataPoint model
</li>
</ul>
</html>"));
end PressurizationData;
