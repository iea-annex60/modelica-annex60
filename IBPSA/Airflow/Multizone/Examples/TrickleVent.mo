within IBPSA.Airflow.Multizone.Examples;
model TrickleVent
  "Model with a trickle vent modelled using the TableData models"
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Air;

  BoundaryConditions.WeatherData.ReaderTMY3       weaDat(filNam=
        Modelica.Utilities.Files.loadResource(
        "modelica://IBPSA/Resources/weatherdata/USA_CA_San.Francisco.Intl.AP.724940_TMY3.mos"))
    annotation (Placement(transformation(extent={{-88,-38},{-68,-18}})));
  Fluid.Sources.Outside_CpLowRise       west(
    redeclare package Medium = Medium,
    s=5,
    azi=IBPSA.Types.Azimuth.W,
    Cp0=0.6,
    nPorts=1) "Model with outside conditions"
    annotation (Placement(transformation(extent={{86,-38},{66,-18}})));
  Fluid.Sources.Outside_CpLowRise east(
    redeclare package Medium = Medium,
    s=5,
    azi=IBPSA.Types.Azimuth.E,
    Cp0=0.6,
    nPorts=1) "Model with outside conditions"
    annotation (Placement(transformation(extent={{-56,-38},{-36,-18}})));
  Fluid.MixingVolumes.MixingVolume       Room(
    redeclare package Medium = Medium,
    V=2.5*5*5,
    nPorts=2,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    m_flow_nominal=0.01)
    annotation (Placement(transformation(extent={{4,-22},{24,-2}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow preHea
    "Prescribed heat flow" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={6,52})));
  Modelica.Blocks.Continuous.LimPID con(
    Td=10,
    yMax=1,
    yMin=-1,
    Ti=60,
    controllerType=Modelica.Blocks.Types.SimpleController.P,
    k=5) "Controller to maintain volume temperature"
    annotation (Placement(transformation(extent={{-54,62},{-34,82}})));
  Modelica.Blocks.Sources.Constant TSet(k=293.15) "Temperature set point"
    annotation (Placement(transformation(extent={{-84,62},{-64,82}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temSen
    "Temperature sensor" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        origin={-76,42})));
  Modelica.Blocks.Math.Gain gain(k=3000)
    annotation (Placement(transformation(extent={{-22,62},{-2,82}})));
  TableData_m_flow tabdat_M(redeclare package Medium = Medium, table=[-50,-0.08709;
        -25,-0.06158; -10,-0.03895; -5,-0.02754; -3,-0.02133; -2,-0.01742; -1,-0.01232;
        0,0; 1,0.01232; 2,0.01742; 3,0.02133; 4.5,0.02613; 50,0.02614])
    "Self regulating trickle vent"
    annotation (Placement(transformation(extent={{-24,-38},{-4,-18}})));
  TableData_V_flow tabdat_V(redeclare package Medium = Medium, table=[-50,-0.104508;
        -25,-0.073896; -10,-0.04674; -5,-0.033048; -3,-0.025596; -2,-0.020904;
        -1,-0.014784; 0,0; 1,0.014784; 2,0.020904; 3,0.025596; 4.5,0.031356; 50,
        0.031368]) "Self regulating trickle vent"
    annotation (Placement(transformation(extent={{36,-38},{56,-18}})));
equation
  connect(weaDat.weaBus, west.weaBus) annotation (Line(
      points={{-68,-28},{-60,-28},{-60,-62},{92,-62},{92,-27.8},{86,-27.8}},
      color={255,204,51},
      thickness=0.5));
  connect(east.weaBus, weaDat.weaBus) annotation (Line(
      points={{-56,-27.8},{-56,-28},{-68,-28}},
      color={255,204,51},
      thickness=0.5));
  connect(TSet.y,con. u_s) annotation (Line(
      points={{-63,72},{-56,72}},
      color={0,0,127}));
  connect(temSen.T,con. u_m) annotation (Line(
      points={{-66,42},{-44,42},{-44,60}},
      color={0,0,127}));
  connect(gain.u,con. y) annotation (Line(
      points={{-24,72},{-33,72}},
      color={0,0,127}));
  connect(gain.y,preHea. Q_flow) annotation (Line(
      points={{-1,72},{6,72},{6,62}},
      color={0,0,127}));
  connect(Room.heatPort, temSen.port) annotation (Line(points={{4,-12},{-6,-12},
          {-6,24},{-86,24},{-86,42}}, color={191,0,0}));
  connect(preHea.port, Room.heatPort) annotation (Line(points={{6,42},{6,24},{
          -2,24},{-2,-12},{4,-12}}, color={191,0,0}));
  connect(east.ports[1], tabdat_M.port_a)
    annotation (Line(points={{-36,-28},{-24,-28}}, color={0,127,255}));
  connect(tabdat_M.port_b, Room.ports[1])
    annotation (Line(points={{-4,-28},{12,-28},{12,-22}}, color={0,127,255}));
  connect(tabdat_V.port_a, Room.ports[2]) annotation (Line(points={{36,-28},{14,
          -28},{14,-22},{16,-22}}, color={0,127,255}));
  connect(tabdat_V.port_b, west.ports[1])
    annotation (Line(points={{56,-28},{66,-28}}, color={0,127,255}));
  annotation (__Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/Examples/OneEffectiveAirLeakageArea.mos"
        "Simulate and plot"),
        experiment(
      StopTime=2592000,
      Interval=600,
      Tolerance=1e-08,
      __Dymola_Algorithm="Dassl"),
    Documentation(info="<html>
<p>
This model illustrates the use of the TableData models of the airlfow.multizone package to model self regulating inlet vents. 
The models are connected to a common volume/room on one side and to outside conditions on the other side (east and west oriëntation respectively).

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
