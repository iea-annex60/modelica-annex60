within IBPSA.Airflow.Multizone.Examples;
model TrickleVent
  "Model with a trickle vent modelled using the TableData models"
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Air;

  BoundaryConditions.WeatherData.ReaderTMY3       weaDat(filNam=
        Modelica.Utilities.Files.loadResource(
        "modelica://IBPSA/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos"))
    "Weather data reader"
    annotation (Placement(transformation(extent={{-100,-40},{-80,-20}})));
  Fluid.Sources.Outside_CpLowRise       west(
    redeclare package Medium = Medium,
    s=5,
    azi=IBPSA.Types.Azimuth.W,
    Cp0=0.6,
    nPorts=1) "Model with outside conditions"
    annotation (Placement(transformation(extent={{80,-40},{60,-20}})));
  Fluid.Sources.Outside_CpLowRise east(
    redeclare package Medium = Medium,
    s=5,
    azi=IBPSA.Types.Azimuth.E,
    Cp0=0.6,
    nPorts=1) "Model with outside conditions"
    annotation (Placement(transformation(extent={{-60,-40},{-40,-20}})));
  Fluid.MixingVolumes.MixingVolume       room(
    redeclare package Medium = Medium,
    V=2.5*5*5,
    nPorts=2,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    m_flow_nominal=0.01) "Room model"
    annotation (Placement(transformation(extent={{0,-20},{20,0}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow preHea
    "Prescribed heat flow" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={10,50})));
  Modelica.Blocks.Continuous.LimPID con(
    Td=10,
    yMax=1,
    yMin=-1,
    Ti=60,
    controllerType=Modelica.Blocks.Types.SimpleController.P,
    k=5) "Controller to maintain volume temperature"
    annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
  Modelica.Blocks.Sources.Constant TSet(k=293.15) "Temperature set point"
    annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temSen
    "Temperature sensor" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        origin={-70,40})));
  Modelica.Blocks.Math.Gain gain(k=3000) "Gain block"
    annotation (Placement(transformation(extent={{-20,60},{0,80}})));
  TableData_m_flow tabdat_M(redeclare package Medium = Medium, table=[-50,-0.08709;
        -25,-0.06158; -10,-0.03895; -5,-0.02754; -3,-0.02133; -2,-0.01742; -1,-0.01232;
        0,0; 1,0.01232; 2,0.01742; 3,0.02133; 4.5,0.02613; 50,0.02614])
    "Self regulating trickle vent"
    annotation (Placement(transformation(extent={{-30,-40},{-10,-20}})));
  TableData_V_flow tabdat_V(redeclare package Medium = Medium, table=[-50,-0.104508;
        -25,-0.073896; -10,-0.04674; -5,-0.033048; -3,-0.025596; -2,-0.020904;
        -1,-0.014784; 0,0; 1,0.014784; 2,0.020904; 3,0.025596; 4.5,0.031356; 50,
        0.031368]) "Self regulating trickle vent"
    annotation (Placement(transformation(extent={{30,-40},{50,-20}})));
equation
  connect(weaDat.weaBus, west.weaBus) annotation (Line(
      points={{-80,-30},{-60,-30},{-60,-60},{80,-60},{80,-29.8}},
      color={255,204,51},
      thickness=0.5));
  connect(east.weaBus, weaDat.weaBus) annotation (Line(
      points={{-60,-29.8},{-60,-30},{-80,-30}},
      color={255,204,51},
      thickness=0.5));
  connect(TSet.y,con. u_s) annotation (Line(
      points={{-79,70},{-62,70}},
      color={0,0,127}));
  connect(temSen.T,con. u_m) annotation (Line(
      points={{-60,40},{-50,40},{-50,58}},
      color={0,0,127}));
  connect(gain.u,con. y) annotation (Line(
      points={{-22,70},{-39,70}},
      color={0,0,127}));
  connect(gain.y,preHea. Q_flow) annotation (Line(
      points={{1,70},{10,70},{10,60}},
      color={0,0,127}));
  connect(room.heatPort, temSen.port) annotation (Line(points={{0,-10},{0,20},{
          -80,20},{-80,40}},          color={191,0,0}));
  connect(preHea.port,room. heatPort) annotation (Line(points={{10,40},{10,20},
          {0,20},{0,-10}},          color={191,0,0}));
  connect(east.ports[1], tabdat_M.port_a)
    annotation (Line(points={{-40,-30},{-30,-30}}, color={0,127,255}));
  connect(tabdat_M.port_b,room. ports[1])
    annotation (Line(points={{-10,-30},{8,-30},{8,-20}},  color={0,127,255}));
  connect(tabdat_V.port_a,room. ports[2]) annotation (Line(points={{30,-30},{14,
          -30},{14,-20},{12,-20}}, color={0,127,255}));
  connect(tabdat_V.port_b, west.ports[1])
    annotation (Line(points={{50,-30},{60,-30}}, color={0,127,255}));
  annotation (__Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/Examples/OneEffectiveAirLeakageArea.mos"
        "Simulate and plot"),
        experiment(
      StopTime=2592000,
      Interval=600,
      Tolerance=1e-08),
    Documentation(info="<html>
<p>
This model illustrates the use of the TableData models of the 
Airflow.MultiZone package to model self regulating inlet vents. 
The models are connected to a common volume/room on one side and 
to outside conditions on the other side (east and west oriëntation respectively).
</p>
</html>", revisions="<html>
<ul>
<li>
May 03, 2021 by Klaas De Jonge:<br/>
Added example for simulating a trickle vent using the TableData models
</li>
</ul>
</html>"));
end TrickleVent;
