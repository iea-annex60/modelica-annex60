within IBPSA.Airflow.Multizone.Examples;
model Powerlaw "Model with powerlaw models"
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Air;
  Powerlaw_m_flow                 powlaw_M(
    redeclare package Medium = Medium,
    m=0.59,
    C=3.33e-5) "Mass flow rate based on powerlaw, direct input for m and C"
             annotation (Placement(transformation(extent={{-2,74},{18,94}})));
  IBPSA.Fluid.Sources.Boundary_pT roo1(
    redeclare package Medium = Medium,
    use_p_in=true,
    nPorts=1,
    T=278.15) annotation (Placement(transformation(extent={{-60,74},{-40,94}})));
  IBPSA.Fluid.Sources.Boundary_pT roo2(
    redeclare package Medium = Medium,
    use_p_in=true,
    nPorts=1,
    T=293.15) annotation (Placement(transformation(
        origin={68,84},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  Modelica.Blocks.Sources.Ramp Ramp1(
    duration=0.5,
    height=2,
    offset=-1,
    startTime=0.25) annotation (Placement(transformation(extent={{2,-80},{22,
            -60}})));
  Modelica.Blocks.Sources.Constant Pre(k=100000) annotation (Placement(
        transformation(extent={{-100,-60},{-80,-40}})));
  Modelica.Blocks.Math.Add Add1 annotation (Placement(transformation(extent={{40,-60},
            {60,-40}})));
  IBPSA.Fluid.Sensors.DensityTwoPort den1(
    redeclare package Medium = Medium,
    m_flow_nominal=0.1,
    initType=Modelica.Blocks.Types.Init.InitialState) "Density sensor"
    annotation (Placement(transformation(extent={{-32,74},{-12,94}})));
  IBPSA.Fluid.Sensors.DensityTwoPort den2(
    redeclare package Medium = Medium,
    m_flow_nominal=0.1,
    initType=Modelica.Blocks.Types.Init.InitialState) "Density sensor"
    annotation (Placement(transformation(extent={{28,74},{48,94}})));
  Powerlaw_V_flow powlaw_V(
    redeclare package Medium = Medium,
    m=0.59,
    C=3.33e-5/1.2)
    "Volume flow rate based on powerlaw, direct input for m and C"
    annotation (Placement(transformation(extent={{-2,42},{18,62}})));
  Powerlaw_1Datapoint powlaw_1dat(
    redeclare package Medium = Medium,
    m=0.59,
    dP1(displayUnit="Pa") = 50,
    m1_flow=1.2/3600)
    "Mass flow rate based on powerlaw, input of m and 1 test data point."
    annotation (Placement(transformation(extent={{-2,10},{18,30}})));
  Powerlaw_2Datapoints powlaw_2dat(
    redeclare package Medium = Medium,
    dP1(displayUnit="Pa") = 1,
    m1_flow=0.12/3600,
    dP2(displayUnit="Pa") = 50,
    m2_flow=1.2/3600)
    "Mass flow rate based on powerlaw, input of 2 test data points."
    annotation (Placement(transformation(extent={{-2,-22},{18,-2}})));
equation
  connect(Pre.y, Add1.u1) annotation (Line(points={{-79,-50},{26,-50},{26,-44},
          {38,-44}}, color={0,0,255}));
  connect(Ramp1.y, Add1.u2) annotation (Line(points={{23,-70},{26,-70},{26,-56},
          {38,-56}}, color={0,0,255}));
  connect(Pre.y, roo1.p_in) annotation (Line(
      points={{-79,-50},{-70,-50},{-70,92},{-62,92}},
      color={0,0,127}));
  connect(Add1.y, roo2.p_in) annotation (Line(
      points={{61,-50},{90,-50},{90,76},{80,76}},
      color={0,0,127}));
  connect(roo1.ports[1], den1.port_a) annotation (Line(
      points={{-40,84},{-32,84}},
      color={0,127,255}));
  connect(den1.port_b, powlaw_M.port_a)
    annotation (Line(points={{-12,84},{-2,84}}, color={0,127,255}));
  connect(powlaw_M.port_b, den2.port_a)
    annotation (Line(points={{18,84},{28,84}}, color={0,127,255}));
  connect(den2.port_b, roo2.ports[1]) annotation (Line(
      points={{48,84},{58,84}},
      color={0,127,255}));
  connect(den1.port_b, powlaw_V.port_a) annotation (Line(points={{-12,84},{-8,
          84},{-8,52},{-2,52}}, color={0,127,255}));
  connect(powlaw_V.port_b, den2.port_a) annotation (Line(points={{18,52},{24,52},
          {24,84},{28,84}}, color={0,127,255}));
  connect(den1.port_b, powlaw_1dat.port_a) annotation (Line(points={{-12,84},{
          -10,84},{-10,20},{-2,20}}, color={0,127,255}));
  connect(den1.port_b, powlaw_2dat.port_a)
    annotation (Line(points={{-12,84},{-12,-12},{-2,-12}}, color={0,127,255}));
  connect(powlaw_1dat.port_b, den2.port_a) annotation (Line(points={{18,20},{26,
          20},{26,84},{28,84}}, color={0,127,255}));
  connect(powlaw_2dat.port_b, den2.port_a)
    annotation (Line(points={{18,-12},{28,-12},{28,84}}, color={0,127,255}));
  annotation (
experiment(Tolerance=1e-06, StopTime=1),
    __Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/Examples/Orifice.mos"
        "Simulate and plot"),
    Documentation(info="<html>
<p>
This model demonstrates the use of the 4 Powerlaw models present in the multizone package. The input data is fit so that all models have equivalent output.

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
</html>"));
end Powerlaw;
