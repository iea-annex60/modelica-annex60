within IBPSA.Fluid.HeatExchangers.Validation;
model WetCoilEffectivenessNTU
  "Validation of WetCoilEffectivesnessNTU based on fan coil unit performance data"
  extends Modelica.Icons.Example;
  package MediumAir = IBPSA.Media.Air;
  package MediumWater = IBPSA.Media.Water;
  parameter Modelica.SIunits.SpecificHeatCapacity cpWat = MediumWater.cp_const;
  IBPSA.Fluid.HeatExchangers.WetCoilEffectivesnessNTU cooCoi16(
    redeclare package Medium1 = MediumAir,
    redeclare package Medium2 = MediumWater,
    Q_flow_nominal=1176,
    dp1_nominal=0,
    dp2_nominal=0,
    m1_flow_nominal=485*1.2/3600,
    configuration=IBPSA.Fluid.Types.HeatExchangerConfiguration.CrossFlowUnmixed,
    m2_flow_nominal=1176/2/4180,
    X_w1_nominal=0.010203,
    T_a1_nominal=300.15,
    T_a2_nominal=289.15,
    X_w2_nominal=0)
    "Cooling coil with nominal conditions for cooling at 16 degrees water inlet temperature"
    annotation (Placement(transformation(extent={{10,40},{-10,20}})));

  Sources.MassFlowSource_T bouAirCoo(
    redeclare package Medium = MediumAir,
    use_m_flow_in=true,
    T=273.15 + 27,
    nPorts=1,
    X={0.010203,1 - 0.010203})
    "Air boundary: 27/19 dry/wet bulb temperature"
    annotation (Placement(transformation(extent={{58,14},{38,34}})));
    // fixme: reformat the table so it is readable.
  Modelica.Blocks.Sources.CombiTimeTable cooData(smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
      table=[0,485,1176,2176,2934; 10000,420,1053,1948,2657; 20000,330,865,1600,
        2209; 30000,236,626,1158,1618; 40000,123,389,720,1018; 50000,0,0.001,0.001,
        0.001])
    "Manufacturers data for cooling mode - 1. air flow rate 2.power at 16 degC 3.sensible power at 7 degC 4.total power at 7degC"
    annotation (Placement(transformation(extent={{-80,70},{-60,90}})));
  Sources.MassFlowSource_T bouAirCoo1(
    redeclare package Medium = MediumAir,
    use_m_flow_in=true,
    T=273.15 + 27,
    nPorts=1,
    X={0.010203,1 - 0.010203})
    "Air boundary: 27/19 dry/wet bulb temperature"
    annotation (Placement(transformation(extent={{58,-26},{38,-6}})));
  Modelica.Blocks.Math.Gain gainFloFcu16(k=1/4180/2)
    "Conversion factor from power to kg/s water "
    annotation (Placement(transformation(extent={{-66,38},{-54,50}})));
  Modelica.Blocks.Math.Gain gain(k=1.2/3600)
    "Conversion factor from l/h to kg/s"
    annotation (Placement(transformation(extent={{-2,50},{18,70}})));
  Sources.MassFlowSource_T bouWatCoo16(
    use_m_flow_in=true,
    redeclare package Medium = MediumWater,
    T=273.15 + 16,
    nPorts=1) "Water boundary at 16 degrees"
    annotation (Placement(transformation(extent={{-40,26},{-20,46}})));
  Sources.MassFlowSource_T bouWatCoo7(
    use_m_flow_in=true,
    redeclare package Medium = MediumWater,
    T=273.15 + 7,
    nPorts=1) "Water boundary at 7 degrees"
    annotation (Placement(transformation(extent={{-40,-14},{-20,6}})));
  Modelica.Blocks.Math.Gain gainFloFcu7(k=1/4180/5)
    "Conversion factor from power to kg/s water "
    annotation (Placement(transformation(extent={{-66,-2},{-54,10}})));
  IBPSA.Fluid.HeatExchangers.WetCoilEffectivesnessNTU cooCoi7(
    redeclare package Medium1 = MediumAir,
    redeclare package Medium2 = MediumWater,
    dp1_nominal=0,
    dp2_nominal=0,
    m1_flow_nominal=485*1.2/3600,
    configuration=IBPSA.Fluid.Types.HeatExchangerConfiguration.CrossFlowUnmixed,
    m2_flow_nominal=2934/5/4180,
    Q_flow_nominal=2934,
    X_w1_nominal=0.010203,
    T_a1_nominal=300.15,
    T_a2_nominal=280.15,
    X_w2_nominal=0)
    "Cooling coil with nominal conditions for cooling at 7 degrees water inlet temperature"
    annotation (Placement(transformation(extent={{10,0},{-10,-20}})));

  Sources.Boundary_pT sinAir(redeclare package Medium = MediumAir, nPorts=4)
    "Air sink"
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Sources.Boundary_pT sinWat(redeclare package Medium = MediumWater, nPorts=4)
    "Water sink"
    annotation (Placement(transformation(extent={{100,0},{80,20}})));
  IBPSA.Fluid.HeatExchangers.WetCoilEffectivesnessNTU cooCoi7Param16(
    redeclare package Medium1 = MediumAir,
    redeclare package Medium2 = MediumWater,
    Q_flow_nominal=1176,
    dp1_nominal=0,
    dp2_nominal=0,
    m1_flow_nominal=485*1.2/3600,
    configuration=IBPSA.Fluid.Types.HeatExchangerConfiguration.CrossFlowUnmixed,
    m2_flow_nominal=1176/2/4180,
    X_w1_nominal=0.010203,
    T_a1_nominal=300.15,
    T_a2_nominal=289.15,
    X_w2_nominal=0)
    "Cooling coil with nominal conditions for cooling at 7 degrees water inlet temperature and boundary conditions of 16 degrees water inlet"
    annotation (Placement(transformation(extent={{10,-40},{-10,-60}})));

  Sources.MassFlowSource_T bouWatCoo7v2(
    use_m_flow_in=true,
    redeclare package Medium = MediumWater,
    T=273.15 + 7,
    nPorts=1) "Water boundary at 7 degrees"
    annotation (Placement(transformation(extent={{-40,-54},{-20,-34}})));
  Sources.MassFlowSource_T bouAirCoo2(
    redeclare package Medium = MediumAir,
    use_m_flow_in=true,
    T=273.15 + 27,
    nPorts=1,
    X={0.010203,1 - 0.010203})
    "Air boundary: 27/19 dry/wet bulb temperature"
    annotation (Placement(transformation(extent={{58,-66},{38,-46}})));
  Sources.MassFlowSource_T bouWatCoo7v3(
    use_m_flow_in=true,
    redeclare package Medium = MediumWater,
    nPorts=1,
    T=273.15 + 7) "Water boundary at 16 degrees"
    annotation (Placement(transformation(extent={{-40,-88},{-20,-68}})));
  Sources.MassFlowSource_T bouAirCoo3(
    redeclare package Medium = MediumAir,
    use_m_flow_in=true,
    T=273.15 + 27,
    X={0.010203,1 - 0.010203},
    nPorts=1)
    "Air boundary: 27/19 dry/wet bulb temperature"
    annotation (Placement(transformation(extent={{58,-100},{38,-80}})));
  WetCoilEffectivesnessNTU cooCoi7Reverse(
    dp1_nominal=0,
    dp2_nominal=0,
    configuration=IBPSA.Fluid.Types.HeatExchangerConfiguration.CrossFlowUnmixed,
    Q_flow_nominal=2934,
    X_w1_nominal=0,
    X_w2_nominal=0.010203,
    redeclare package Medium1 = MediumWater,
    redeclare package Medium2 = MediumAir,
    m1_flow_nominal=2934/5/4180,
    m2_flow_nominal=485*1.2/3600,
    T_a1_nominal=280.15,
    T_a2_nominal=300.15)
    "Reverse configuration for verifying symmetry of the implementation"
    annotation (Placement(transformation(extent={{-10,-94},{10,-74}})));

equation
  // fixme: enable or delete the assertion
  //assert(abs(cooCoi7.Q1_flow-cooCoi7Reverse.Q1_flow)+
  //       abs(cooCoi7.mWat1_flow-cooCoi7Reverse.mWat2_flow)<1e-10,
  //      "The models cooCoi7 and cooCoi7Reverse do
  //      not return identical results as required.");
  connect(gainFloFcu16.u, cooData.y[2])
    annotation (Line(points={{-67.2,44},{-74,44},{-74,60},{-50,60},{-50,80},{-59,
          80}},                                      color={0,0,127}));
  connect(gain.u, cooData.y[1])
    annotation (Line(points={{-4,60},{-50,60},{-50,80},{-59,80}},
                                                  color={0,0,127}));
  connect(gain.y, bouAirCoo.m_flow_in)
    annotation (Line(points={{19,60},{68,60},{68,32},{60,32}},
                                                          color={0,0,127}));
  connect(gainFloFcu16.y, bouWatCoo16.m_flow_in) annotation (Line(points={{-53.4,
          44},{-42,44}},             color={0,0,127}));
  connect(gainFloFcu7.y,bouWatCoo7. m_flow_in)
    annotation (Line(points={{-53.4,4},{-42,4}},     color={0,0,127}));
  connect(cooCoi16.port_a2, bouWatCoo16.ports[1])
    annotation (Line(points={{-10,36},{-20,36}},  color={0,127,255}));
  connect(cooCoi16.port_b2, sinWat.ports[1])
    annotation (Line(points={{10,36},{80,36},{80,13}},  color={0,127,255}));
  connect(cooCoi16.port_a1, bouAirCoo.ports[1])
    annotation (Line(points={{10,24},{38,24}},   color={0,127,255}));
  connect(cooCoi16.port_b1, sinAir.ports[1]) annotation (Line(points={{-10,24},{
          -80,24},{-80,13}},   color={0,127,255}));
  connect(cooCoi7.port_a2, bouWatCoo7.ports[1])
    annotation (Line(points={{-10,-4},{-20,-4}},  color={0,127,255}));
  connect(cooCoi7.port_b2, sinWat.ports[2])
    annotation (Line(points={{10,-4},{80,-4},{80,11}},   color={0,127,255}));
  connect(bouAirCoo1.ports[1], cooCoi7.port_a1)
    annotation (Line(points={{38,-16},{10,-16}}, color={0,127,255}));
  connect(cooCoi7.port_b1, sinAir.ports[2]) annotation (Line(points={{-10,-16},{
          -80,-16},{-80,11}},
                           color={0,127,255}));
  connect(gainFloFcu7.u, cooData.y[4]) annotation (Line(points={{-67.2,4},{-74,4},
          {-74,60},{-50,60},{-50,80},{-59,80}},
                                     color={0,0,127}));
  connect(bouWatCoo7v2.ports[1], cooCoi7Param16.port_a2)
    annotation (Line(points={{-20,-44},{-10,-44}}, color={0,127,255}));
  connect(bouWatCoo7v2.m_flow_in, gainFloFcu7.y) annotation (Line(points={{-42,-36},
          {-50,-36},{-50,4},{-53.4,4}},color={0,0,127}));
  connect(bouAirCoo2.ports[1], cooCoi7Param16.port_a1)
    annotation (Line(points={{38,-56},{10,-56}}, color={0,127,255}));
  connect(sinWat.ports[3], cooCoi7Param16.port_b2) annotation (Line(points={{80,9},{
          80,-44},{10,-44}},           color={0,127,255}));
  connect(cooCoi7Param16.port_b1, sinAir.ports[3]) annotation (Line(points={{-10,-56},
          {-80,-56},{-80,9}},                color={0,127,255}));
  connect(cooCoi7Reverse.port_a1, bouWatCoo7v3.ports[1])
    annotation (Line(points={{-10,-78},{-20,-78}}, color={0,127,255}));
  connect(cooCoi7Reverse.port_a2, bouAirCoo3.ports[1])
    annotation (Line(points={{10,-90},{38,-90}}, color={0,127,255}));
  connect(cooCoi7Reverse.port_b2, sinAir.ports[4])
    annotation (Line(points={{-10,-90},{-80,-90},{-80,7}}, color={0,127,255}));
  connect(cooCoi7Reverse.port_b1, sinWat.ports[4])
    annotation (Line(points={{10,-78},{80,-78},{80,7}}, color={0,127,255}));
  connect(bouWatCoo7v3.m_flow_in, gainFloFcu7.y) annotation (Line(points={{-42,-70},
          {-50,-70},{-50,4},{-53.4,4}},
                                  color={0,0,127}));
  connect(gain.y, bouAirCoo1.m_flow_in) annotation (Line(points={{19,60},{68,60},
          {68,-8},{60,-8}}, color={0,0,127}));
  connect(gain.y, bouAirCoo2.m_flow_in) annotation (Line(points={{19,60},{68,60},
          {68,-48},{60,-48}}, color={0,0,127}));
  connect(gain.y, bouAirCoo3.m_flow_in) annotation (Line(points={{19,60},{68,60},
          {68,-82},{60,-82}}, color={0,0,127}));
  annotation (experiment(
      StopTime=50000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Lsodar"), __Dymola_Commands(file="Resources/Scripts/Dymola/Fluid/HeatExchangers/Validation/WetCoilEffectivenessNTU.mos"
        "Simulate and plot"),
    Documentation(info="<html>
<p>
fixme: complete info section.
</p>
</html>", revisions="<html>
<p>
fixme: add revision section
</p>
</html>"));
end WetCoilEffectivenessNTU;
