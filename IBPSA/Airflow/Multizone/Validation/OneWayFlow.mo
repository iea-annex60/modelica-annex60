within IBPSA.Airflow.Multizone.Validation;
model OneWayFlow
  "Validation model to verify one way flow implementation"
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Specialized.Air.PerfectGas;


protected
  parameter Integer nTested=8 "Number of tested flow elements";

//Test Data
//Headers: ["dP","ELA_FlowRate","ORI_FlowRate","1DatPoint_FlowRate","2DatPoint_FlowRate","PowLaw_M_FlowRate","PowLaw_V_FlowRate","TabDat_M","TabDat_V_FlowRate"]
  parameter Real TestData[:,nTested+1]=[
      -50,-0.0838,-0.0658,-0.0672,-0.0609,-0.0707,-0.0851,-0.0871,-0.105;
      -40,-0.0725,-0.0589,-0.0601,-0.055,-0.0632,-0.0762,-0.0769,-0.0926;
      -25,-0.0534,-0.0466,-0.0475,-0.0443,-0.05,-0.0602,-0.0616,-0.0741;
      -10,-0.0294,-0.0294,-0.03,-0.029,-0.0316,-0.0381,-0.039,-0.0469;
      -5,-0.0188,-0.0208,-0.0212,-0.0211,-0.0224,-0.0269,-0.0275,-0.0332;
      -1,-0.00659,-0.00931,-0.0095,-0.01,-0.01,-0.012,-0.0123,-0.0148;
      0,0,0,0,0,0,0,0,0;
      1,0.00659,0.00931,0.0095,0.01,0.01,0.012,0.0123,0.0148;
      5,0.0188,0.0208,0.0212,0.0211,0.0224,0.0269,0.0261,0.0315;
      10,0.0294,0.0294,0.03,0.029,0.0316,0.0381,0.0261,0.0315;
      25,0.0534,0.0466,0.0475,0.0443,0.05,0.0602,0.0261,0.0315;
      40,0.0725,0.0589,0.0601,0.055,0.0632,0.0762,0.0261,0.0315;
      50,0.0838,0.0658,0.0672,0.0609,0.0707,0.0851,0.0261,0.0315]
      "Steady state CONTAM results with specific pressure difference for similar flow models";

  //Boundary conditions
  Fluid.Sources.Boundary_pT bouA(
    redeclare package Medium = Medium,
    use_p_in=true,
    T=293.15,
    nPorts=8) annotation (Placement(transformation(extent={{-100,-20},{-80,0}})));
  Fluid.Sources.Boundary_pT bouB(
    redeclare package Medium = Medium,
    use_p_in=true,
    T=293.15,
    nPorts=8) annotation (Placement(transformation(extent={{120,-20},{100,0}})));
  Modelica.Blocks.Sources.Ramp ramp_min50_50pa(
    duration=500,
    height=100,
    offset=-50) "Block that generates a ramp signal from -50 to +50"
     annotation (Placement(transformation(extent={{-140,20},{-120,40}})));
  Modelica.Blocks.Sources.Constant AmbP(k=101325)
  "Assumed ambient pressure" annotation (Placement(transformation(extent={{160,20},{140,40}})));
  Modelica.Blocks.Math.Sum sum(nin=2) annotation (Placement(transformation(
        extent={{-2,-2},{2,2}},
        rotation=0,
        origin={-108,-2})));

  //Flow models
  EffectiveAirLeakageArea ela(
    redeclare package Medium = Medium,
    dpRat=10,
    CDRat=0.6,
    L=0.01)
    "EffectiveAirLeakageArea" annotation (Placement(transformation(extent={{-40,120},
            {-20,140}})));
  Orifice ori(
    redeclare package Medium = Medium,
    A=0.01,
    CD=0.6) "Orifice" annotation (Placement(transformation(extent={{-40,80},{-20,100}})));
  Powerlaw_1Datapoint powlaw_1dat(
    redeclare package Medium = Medium,
    dP1(displayUnit="Pa") = 4,
    m1_flow=0.019)
    "Powerlaw_1Datapoint"  annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  Powerlaw_2Datapoints powlaw_2dat(
    redeclare package Medium = Medium,
    dP1(displayUnit="Pa") = 4,
    m1_flow=0.019,
    dP2(displayUnit="Pa") = 10,
    m2_flow=0.029)
    "Powerlaw_2Datapoints" annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Powerlaw_m_flow powlaw_M(
    redeclare package Medium = Medium,
    m=0.5,
    C=0.01)
    "Powerlaw_m_flow" annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
  Powerlaw_V_flow powlaw_V(
    redeclare package Medium = Medium,
    m=0.5,
    C=0.01)
    "Powerlaw_V_flow" annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
  TableData_m_flow tabdat_M(
   redeclare package Medium = Medium,
    table=[-50,-0.08709;
        -25,-0.06158;
        -10,-0.03895;
        -5,-0.02754;
        -3,-0.02133;
        -2,-0.01742;
        -1,-0.01232;
        0,0;
        1,0.01232;
        2,0.01742;
        3,0.02133;
        4.5,0.02613;
        50,0.02614])
    "TableData_m_flow" annotation (Placement(transformation(extent={{-40,-120},{-20,-100}})));
  TableData_V_flow tabdat_V(
    redeclare package Medium = Medium,
    table=[-50,-0.08709;
        -25,-0.06158;
        -10,-0.03895;
        -5,-0.02754;
        -3,-0.02133;
        -2,-0.01742;
        -1,-0.01232;
        0,0;
        1,0.01232;
        2,0.01742;
        3,0.02133;
        4.5,0.02613;
        50,0.02614])
    "TableData_V_flow" annotation (Placement(transformation(extent={{-40,-160},{-20,-140}})));

  //Mass flow sensors
  Fluid.Sensors.MassFlowRate Sen_ela(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,120},{20,140}})));
  Fluid.Sensors.MassFlowRate Sen_ori(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,80},{20,100}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_1dat(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,40},{20,60}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_2dat(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_M(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,-40},{20,-20}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_V(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,-80},{20,-60}})));
  Fluid.Sensors.MassFlowRate Sen_tabdat_M(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,-120},{20,-100}})));
  Fluid.Sensors.MassFlowRate Sen_tabdat_V(redeclare package Medium = Medium)
  "Mass flow rate sensor" annotation (Placement(transformation(extent={{0,-160},{20,-140}})));

  //Checking the data
  Modelica.Blocks.Tables.CombiTable1D IntTestData(
    table=TestData,
    columns=2:9,
    smoothness=Modelica.Blocks.Types.Smoothness.MonotoneContinuousDerivative2)
    "Table with CONTAM simulation results for comparisson" annotation (Placement(transformation(extent={{-80,160},{-60,180}})));
  Modelica.Blocks.Routing.Replicator replicator(nout=IntTestData.n) annotation (Placement(transformation(extent={{-98,164},{-88,176}})));

public
  Modelica.SIunits.PressureDifference dP "Pressure difference over the tested elements";
  Modelica.SIunits.MassFlowRate[nTested]  m_flow_data "simulated mass flow of each flow element";
  Modelica.SIunits.MassFlowRate [nTested] m_flow_testdata "mass flow of each flow element of CONTAM simulation";
equation

  connect(ramp_min50_50pa.y, sum.u[1]) annotation (Line(points={{-119,30},{-116,
          30},{-116,-2},{-110.4,-2},{-110.4,-2.2}},
                                             color={0,0,127}));
  connect(AmbP.y, sum.u[2]) annotation (Line(points={{139,30},{132,30},{132,156},
          {-114,156},{-114,-2},{-110.4,-2},{-110.4,-1.8}},
                                   color={0,0,127}));
  connect(sum.y, bouA.p_in) annotation (Line(points={{-105.8,-2},{-102,-2}},      color={0,0,127}));
  connect(AmbP.y, bouB.p_in) annotation (Line(points={{139,30},{132,30},{132,-2},{122,-2}},
                                                             color={0,0,127}));
  connect(ramp_min50_50pa.y, replicator.u) annotation (Line(points={{-119,30},{-116,30},{-116,170},{-99,170}},
                                                             color={0,0,127}));
  connect(replicator.y, IntTestData.u) annotation (Line(points={{-87.5,170},{-82,170}},
                                                   color={0,0,127}));
  connect(bouA.ports[1],ela. port_a) annotation (Line(points={{-80,-6.5},{-62,-6.5},
          {-62,130},{-40,130}},
                              color={0,127,255}));
  connect(bouA.ports[2], ori.port_a) annotation (Line(points={{-80,-7.5},{-60,-7.5},
          {-60,90},{-40,90}}, color={0,127,255}));
  connect(bouA.ports[3], powlaw_1dat.port_a) annotation (Line(points={{-80,-8.5},
          {-58,-8.5},{-58,50},{-40,50}}, color={0,127,255}));
  connect(bouA.ports[4], powlaw_2dat.port_a) annotation (Line(points={{-80,-9.5},
          {-56,-9.5},{-56,10},{-40,10}}, color={0,127,255}));
  connect(bouA.ports[5], powlaw_M.port_a) annotation (Line(points={{-80,-10.5},{
          -56,-10.5},{-56,-30},{-40,-30}},
                                      color={0,127,255}));
  connect(bouA.ports[6], powlaw_V.port_a) annotation (Line(points={{-80,-11.5},{
          -58,-11.5},{-58,-70},{-40,-70}},
                                      color={0,127,255}));
  connect(bouA.ports[7],tabdat_M. port_a) annotation (Line(points={{-80,-12.5},{
          -60,-12.5},{-60,-110},{-40,-110}},
                                      color={0,127,255}));
  connect(bouA.ports[8],tabdat_V. port_a) annotation (Line(points={{-80,-13.5},{
          -62,-13.5},{-62,-150},{-40,-150}},
                                        color={0,127,255}));
  connect(ela.port_b, Sen_ela.port_a)    annotation (Line(points={{-20,130},{0,130}},
                                              color={0,127,255}));
  connect(Sen_ela.port_b, bouB.ports[1]) annotation (Line(points={{20,130},{60,130},
          {60,-6.5},{100,-6.5}},color={0,127,255}));
  connect(ori.port_b, Sen_ori.port_a)    annotation (Line(points={{-20,90},{0,90}},color={0,127,255}));
  connect(Sen_ori.port_b, bouB.ports[2]) annotation (Line(points={{20,90},{58,90},
          {58,-7.5},{100,-7.5}},color={0,127,255}));
  connect(powlaw_1dat.port_b, Sen_powlaw_1dat.port_a) annotation (Line(points={{-20,50},{0,50}},color={0,127,255}));
  connect(Sen_powlaw_1dat.port_b, bouB.ports[3]) annotation (Line(points={{20,50},
          {56,50},{56,-8.5},{100,-8.5}},color={0,127,255}));
  connect(powlaw_2dat.port_b, Sen_powlaw_2dat.port_a) annotation (Line(points={{-20,10},{0,10}},color={0,127,255}));
  connect(Sen_powlaw_2dat.port_b, bouB.ports[4]) annotation (Line(points={{20,10},
          {54,10},{54,-9.5},{100,-9.5}},color={0,127,255}));
  connect(powlaw_M.port_b, Sen_powlaw_M.port_a) annotation (Line(points={{-20,-30},{0,-30}},color={0,127,255}));
  connect(Sen_powlaw_M.port_b, bouB.ports[5]) annotation (Line(points={{20,-30},
          {54,-30},{54,-10.5},{100,-10.5}},
                                         color={0,127,255}));
  connect(powlaw_V.port_b, Sen_powlaw_V.port_a) annotation (Line(points={{-20,-70},{0,-70}},color={0,127,255}));
  connect(Sen_powlaw_V.port_b, bouB.ports[6]) annotation (Line(points={{20,-70},
          {56,-70},{56,-11.5},{100,-11.5}},
                                color={0,127,255}));
  connect(tabdat_M.port_b, Sen_tabdat_M.port_a) annotation (Line(points={{-20,-110},{0,-110}},
                                                color={0,127,255}));
  connect(Sen_tabdat_M.port_b, bouB.ports[7]) annotation (Line(points={{20,-110},
          {58,-110},{58,-12.5},{100,-12.5}},
                                         color={0,127,255}));
  connect(tabdat_V.port_b, Sen_tabdat_V.port_a) annotation (Line(points={{-20,-150},{0,-150}},  color={0,127,255}));
  connect(Sen_tabdat_V.port_b, bouB.ports[8]) annotation (Line(points={{20,-150},
          {60,-150},{60,-13.5},{100,-13.5}},
                                          color={0,127,255}));

  dP=ela.dp;

  Sen_ela.m_flow= m_flow_data[1];
  Sen_ori.m_flow= m_flow_data[2];
  Sen_powlaw_1dat.m_flow= m_flow_data[3];
  Sen_powlaw_2dat.m_flow= m_flow_data[4];
  Sen_powlaw_M.m_flow= m_flow_data[5];
  Sen_powlaw_V.m_flow= m_flow_data[6];
  Sen_tabdat_M.m_flow= m_flow_data[7];
  Sen_tabdat_V.m_flow= m_flow_data[8];

  IntTestData.y= m_flow_testdata;

  annotation (Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-140,-180},{160,200}})),
    experiment(
      StopTime=500,
      Interval=1,
      Tolerance=1e-06),
      __Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/Validation/OneWayFlow.mos"
        "Simulate and plot"),
    Documentation(revisions="<html>
<ul>
<li>Apr 6, 2021, 2020, by Klaas De Jonge (UGent):<br/>
First implementation</li>
</ul>
</html>
"));
end OneWayFlow;
