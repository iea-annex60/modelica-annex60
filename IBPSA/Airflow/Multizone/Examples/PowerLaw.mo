within IBPSA.Airflow.Multizone.Examples;
model PowerLaw "Model with powerlaw models"
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Air;
  Powerlaw_m_flow                 powlaw_M(
    redeclare package Medium = Medium,
    m=0.59,
    C=3.33e-5) "Mass flow rate based on powerlaw, direct input for m and C"
             annotation (Placement(transformation(extent={{-10,80},{10,100}})));
  IBPSA.Fluid.Sources.Boundary_pT roo1(
    redeclare package Medium = Medium,
    use_p_in=true,
    nPorts=1,
    T=278.15) "Room 1"
              annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
  IBPSA.Fluid.Sources.Boundary_pT roo2(
    redeclare package Medium = Medium,
    use_p_in=true,
    nPorts=1,
    T=293.15) "Room 2"
              annotation (Placement(transformation(
        origin={90,70},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica.Blocks.Sources.Ramp ramp1(
    duration=0.5,
    height=6,
    offset=-3,
    startTime=0.25) "Ramp"
                    annotation (Placement(transformation(extent={{0,-100},{20,-80}})));
  Modelica.Blocks.Sources.Constant pressure(k=100000) "Pressure"
    annotation (Placement(transformation(extent={{-120,-80},{-100,-60}})));
  Modelica.Blocks.Math.Add add "Add"
    annotation (Placement(transformation(extent={{40,-86},{60,-66}})));
  IBPSA.Fluid.Sensors.DensityTwoPort den1(
    redeclare package Medium = Medium,
    m_flow_nominal=0.1,
    tau=0,
    initType=Modelica.Blocks.Types.Init.InitialState) "Density sensor"
    annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
  IBPSA.Fluid.Sensors.DensityTwoPort den2(
    redeclare package Medium = Medium,
    m_flow_nominal=0.1,
    tau=0,
    initType=Modelica.Blocks.Types.Init.InitialState) "Density sensor"
    annotation (Placement(transformation(extent={{40,60},{60,80}})));
  Powerlaw_V_flow powlaw_V(
    redeclare package Medium = Medium,
    m=0.59,
    C=3.33e-5/1.2)
    "Volume flow rate based on powerlaw, direct input for m and C"
    annotation (Placement(transformation(extent={{-10,40},{10,60}})));
  Powerlaw_1Datapoint powlaw_1dat(
    redeclare package Medium = Medium,
    m=0.59,
    dP1(displayUnit="Pa") = 50,
    m1_flow=1.2/3600)
    "Mass flow rate based on powerlaw, input of m and 1 test data point."
    annotation (Placement(transformation(extent={{-10,0},{10,20}})));
  Powerlaw_2Datapoints powlaw_2dat(
    redeclare package Medium = Medium,
    dP1(displayUnit="Pa") = 1,
    m1_flow=0.12/3600,
    dP2(displayUnit="Pa") = 50,
    m2_flow=1.2/3600)
    "Mass flow rate based on powerlaw, input of 2 test data points."
    annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
equation
  connect(pressure.y, add.u1)
    annotation (Line(points={{-99,-70},{38,-70}}, color={0,0,255}));
  connect(ramp1.y, add.u2) annotation (Line(points={{21,-90},{26,-90},{26,-82},{
          38,-82}}, color={0,0,255}));
  connect(pressure.y, roo1.p_in) annotation (Line(points={{-99,-70},{-60,-70},{-60,
          -40},{-120,-40},{-120,78},{-102,78}},
                                        color={0,0,127}));
  connect(add.y, roo2.p_in) annotation (Line(points={{61,-76},{112,-76},{112,62},
          {102,62}}, color={0,0,127}));
  connect(roo1.ports[1], den1.port_a) annotation (Line(
      points={{-80,70},{-60,70}},
      color={0,127,255}));
  connect(den1.port_b, powlaw_M.port_a)
    annotation (Line(points={{-40,70},{-40,90},{-10,90}},
                                                color={0,127,255}));
  connect(powlaw_M.port_b, den2.port_a)
    annotation (Line(points={{10,90},{40,90},{40,70}},
                                               color={0,127,255}));
  connect(den2.port_b, roo2.ports[1]) annotation (Line(
      points={{60,70},{80,70}},
      color={0,127,255}));
  connect(den1.port_b, powlaw_V.port_a) annotation (Line(points={{-40,70},{-40,
          50},{-10,50}},        color={0,127,255}));
  connect(powlaw_V.port_b, den2.port_a) annotation (Line(points={{10,50},{40,50},
          {40,70}},         color={0,127,255}));
  connect(den1.port_b, powlaw_1dat.port_a) annotation (Line(points={{-40,70},{
          -40,10},{-10,10}},         color={0,127,255}));
  connect(den1.port_b, powlaw_2dat.port_a)
    annotation (Line(points={{-40,70},{-40,-30},{-10,-30}},color={0,127,255}));
  connect(powlaw_1dat.port_b, den2.port_a) annotation (Line(points={{10,10},{40,
          10},{40,70}},         color={0,127,255}));
  connect(powlaw_2dat.port_b, den2.port_a)
    annotation (Line(points={{10,-30},{40,-30},{40,70}}, color={0,127,255}));
  annotation (
experiment(
      StopTime=1,
      Tolerance=1e-06),
    __Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/Examples/PowerLaw.mos"
        "Simulate and plot"),
    Documentation(info="<html>
<p>
This model demonstrates the use of the 4 PowerLaw models present in the multizone package. 
The input data is fit so that all models have equivalent output.

The pressure difference across the models changes
between <i>-1</i> Pascal and <i>+1</i> Pascal, which
causes air to flow through the orifice.
</p>
</html>", revisions="<html>
<ul>
<li>
May 3, 2021 by Klaas De Jonge:<br/>
Added example for the four 'Powerlaw' models in the Multizone package.
</li>

</ul>
</html>"),
    Diagram(coordinateSystem(extent={{-120,-120},{120,100}})));
end PowerLaw;
