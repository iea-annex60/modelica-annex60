within IBPSA.Airflow.Multizone.Validation;
model OneWayFlow
  "Validation model to verify one way flow implementation"
  extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Specialized.Air.PerfectGas;

//Test Data

  //Headers: dP,ELA_FlowRate,ORI_FlowRate,PowLaw_M_FlowRate,PowLaw_V_FlowRate,TabDat_M,TabDat_V_FlowRate,2DatPoint_FlowRate,1DatPoint_FlowRate
  //parameter String Headers[:,:]=["dP","ELA_FlowRate","ORI_FlowRate","PowLaw_M_FlowRate","PowLaw_V_FlowRate","TabDat_M","TabDat_V_FlowRate","2DatPoint_FlowRate","1DatPoint_FlowRate"];
protected
  parameter Integer nTested=8 "Number of tested flow elements";
  parameter Real TestData[:,:]=[-50,-0.0838,-0.0658,-0.0707,-0.0851,-0.0871,-0.105,
      -0.0609,-0.0672; -40,-0.0725,-0.0589,-0.0632,-0.0762,-0.0769,-0.0926,-0.055,
      -0.0601; -25,-0.0534,-0.0466,-0.05,-0.0602,-0.0616,-0.0741,-0.0443,-0.0475;
      -10,-0.0294,-0.0294,-0.0316,-0.0381,-0.039,-0.0469,-0.029,-0.03; -5,-0.0188,
      -0.0208,-0.0224,-0.0269,-0.0275,-0.0332,-0.0211,-0.0212; -1,-0.00659,-0.00931,
      -0.01,-0.012,-0.0123,-0.0148,-0.01,-0.0095; 0,0,0,0,0,0,0,0,0; 1,0.00659,0.00931,
      0.01,0.012,0.0123,0.0148,0.01,0.0095; 5,0.0188,0.0208,0.0224,0.0269,0.0261,
      0.0315,0.0211,0.0212; 10,0.0294,0.0294,0.0316,0.0381,0.0261,0.0315,0.029,0.03;
      25,0.0534,0.0466,0.05,0.0602,0.0261,0.0315,0.0443,0.0475; 40,0.0725,0.0589,
      0.0632,0.0762,0.0261,0.0315,0.055,0.0601; 50,0.0838,0.0658,0.0707,0.0851,0.0261,
      0.0315,0.0609,0.0672] "CONTAM results of simulations with specific pressure difference for similar flow models";

  //Boundary condition
  Fluid.Sources.Boundary_pT bouA(
    redeclare package Medium = Medium,
    use_p_in=true,
    use_T_in=true,
    nPorts=8) annotation (Placement(transformation(extent={{-80,-14},{-60,6}})));
  Fluid.Sources.Boundary_pT bouB(
    redeclare package Medium = Medium,
    use_p_in=true,
    use_T_in=true,
    nPorts=8) annotation (Placement(transformation(extent={{80,-14},{60,6}})));
  Modelica.Blocks.Sources.Ramp ramp_min50_50pa(
    duration=500,
    height=100,
    offset=-50)
    annotation (Placement(transformation(extent={{-140,20},{-120,40}})));
  Modelica.Blocks.Sources.Constant AmbP(k=101325)
    annotation (Placement(transformation(extent={{120,20},{100,40}})));
  Modelica.Blocks.Sources.Constant Ta(k=20 + 273.15)
    annotation (Placement(transformation(extent={{-140,-10},{-120,10}})));
  Modelica.Blocks.Sources.Constant Tb(k=20 + 273.15)
    annotation (Placement(transformation(extent={{120,-10},{100,10}})));
  Modelica.Blocks.Math.Sum sum(nin=2) annotation (Placement(transformation(
        extent={{-2,-2},{2,2}},
        rotation=270,
        origin={-112,22})));

  //Flow models
  EffectiveAirLeakageArea ela(
    redeclare package Medium = Medium,
    dpRat=10,
    CDRat=0.6,
    L=0.01)
    "EffectiveAirLeakageArea" annotation (Placement(transformation(extent={{-40,96},
            {-20,116}})));

  Orifice ori(
    redeclare package Medium = Medium,
    A=0.01,
    CD=0.6) "Orifice"
    annotation (Placement(transformation(extent={{-40,68},{-20,88}})));
  Powerlaw_1Datapoint powlaw_1dat(
    redeclare package Medium = Medium,
    dP1(displayUnit="Pa") = 4,
    m1_flow=0.019)
               "Powerlaw_1Datapoint"
               annotation (Placement(transformation(extent={{-40,38},{-20,58}})));
  Powerlaw_2Datapoints powlaw_2dat(
    redeclare package Medium = Medium,
    dP1(displayUnit="Pa") = 4,
    m1_flow=0.019,
    dP2(displayUnit="Pa") = 10,
    m2_flow=0.029)
               "Powerlaw_2Datapoints"
               annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Powerlaw_m_flow powlaw_M(
    redeclare package Medium = Medium,
    m=0.5,
    C=0.01) "Powerlaw_m_flow"
            annotation (Placement(transformation(extent={{-40,-28},{-20,-8}})));
  Powerlaw_V_flow powlaw_V(
    redeclare package Medium = Medium,
    m=0.5,
    C=0.01) "Powerlaw_V_flow"
            annotation (Placement(transformation(extent={{-40,-56},{-20,-36}})));
  TableData_m_flow TabDat_M(redeclare package Medium = Medium)
    "TableData_m_flow"
    annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
  TableData_V_flow TabDat_V(redeclare package Medium = Medium)
    "TableData_V_flow"
    annotation (Placement(transformation(extent={{-40,-128},{-20,-108}})));

  //mass flow sensors
  Fluid.Sensors.MassFlowRate Sen_ela(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,96},{12,116}})));
  Fluid.Sensors.MassFlowRate Sen_ori(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,68},{12,88}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_1dat(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,38},{12,58}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_2dat(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,0},{12,20}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_M(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,-28},{12,-8}})));
  Fluid.Sensors.MassFlowRate Sen_powlaw_V(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,-56},{12,-36}})));
  Fluid.Sensors.MassFlowRate Sen_tabdat_M(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,-90},{12,-70}})));
  Fluid.Sensors.MassFlowRate Sen_tabdat_V(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-8,-128},{12,-108}})));

  //Checking the data

                                                         //ndims(TestData)-1

  Modelica.Blocks.Tables.CombiTable1D IntTestData(
    table=TestData,
    columns=2:9,
    smoothness=Modelica.Blocks.Types.Smoothness.MonotoneContinuousDerivative2)
    annotation (Placement(transformation(extent={{-90,66},{-70,86}})));
  Modelica.Blocks.Routing.Replicator replicator(nout=IntTestData.n)
    annotation (Placement(transformation(extent={{-110,70},{-100,82}})));
public
  Modelica.Blocks.Interfaces.RealOutput[nTested]  m_flow_data
    "mass flow of each flow element with order corresponding to test data"
    annotation (Placement(transformation(extent={{116,82},{136,102}})));
  Modelica.Blocks.Interfaces.RealOutput[nTested] m_flow_testdata
    "mass flow of each flow element with order corresponding to test data"
    annotation (Placement(transformation(extent={{116,58},{136,78}})));
equation
  //Boundary condition connections
  connect(bouA.T_in, Ta.y)
    annotation (Line(points={{-82,0},{-119,0}},color={0,0,127}));
  connect(ramp_min50_50pa.y, sum.u[1]) annotation (Line(points={{-119,30},{-112,
          30},{-112,24},{-112.2,24},{-112.2,24.4}},
                                             color={0,0,127}));
  connect(AmbP.y, sum.u[2]) annotation (Line(points={{99,30},{-111.8,30},{
          -111.8,24.4}},           color={0,0,127}));
  connect(sum.y, bouA.p_in)
    annotation (Line(points={{-112,19.8},{-112,4},{-82,4}},
                                                          color={0,0,127}));
  connect(Tb.y, bouB.T_in)
    annotation (Line(points={{99,0},{82,0}}, color={0,0,127}));
  connect(AmbP.y, bouB.p_in)
    annotation (Line(points={{99,30},{92,30},{92,4},{82,4}}, color={0,0,127}));
  //Flow element connections
  connect(bouA.ports[1],ela. port_a) annotation (Line(points={{-60,-0.5},{-54,-0.5},
          {-54,106},{-40,106}},
                              color={0,127,255}));
  connect(bouA.ports[2], ori.port_a) annotation (Line(points={{-60,-1.5},{-50,-1.5},
          {-50,78},{-40,78}}, color={0,127,255}));
  connect(bouA.ports[3], powlaw_1dat.port_a) annotation (Line(points={{-60,-2.5},
          {-50,-2.5},{-50,48},{-40,48}}, color={0,127,255}));
  connect(bouA.ports[4], powlaw_2dat.port_a) annotation (Line(points={{-60,-3.5},
          {-50,-3.5},{-50,10},{-40,10}}, color={0,127,255}));
  connect(bouA.ports[5], powlaw_M.port_a) annotation (Line(points={{-60,-4.5},{-50,
          -4.5},{-50,-18},{-40,-18}}, color={0,127,255}));
  connect(bouA.ports[6], powlaw_V.port_a) annotation (Line(points={{-60,-5.5},{-50,
          -5.5},{-50,-46},{-40,-46}}, color={0,127,255}));
  connect(bouA.ports[7], TabDat_M.port_a) annotation (Line(points={{-60,-6.5},{-50,
          -6.5},{-50,-80},{-40,-80}}, color={0,127,255}));
  connect(bouA.ports[8], TabDat_V.port_a) annotation (Line(points={{-60,-7.5},{-50,
          -7.5},{-50,-118},{-40,-118}}, color={0,127,255}));
  connect(ramp_min50_50pa.y, replicator.u)
    annotation (Line(points={{-119,30},{-111,30},{-111,76}}, color={0,0,127}));
  connect(replicator.y, IntTestData.u)
    annotation (Line(points={{-99.5,76},{-92,76}}, color={0,0,127}));
  connect(ela.port_b, Sen_ela.port_a)
    annotation (Line(points={{-20,106},{-8,106}},
                                              color={0,127,255}));
  connect(Sen_ela.port_b, bouB.ports[1]) annotation (Line(points={{12,106},{28,106},
          {28,-0.5},{60,-0.5}}, color={0,127,255}));
  connect(ori.port_b, Sen_ori.port_a)
    annotation (Line(points={{-20,78},{-8,78}},
                                              color={0,127,255}));
  connect(Sen_ori.port_b, bouB.ports[2]) annotation (Line(points={{12,78},{24,78},
          {24,-1.5},{60,-1.5}}, color={0,127,255}));
  connect(powlaw_1dat.port_b, Sen_powlaw_1dat.port_a)
    annotation (Line(points={{-20,48},{-8,48}},
                                              color={0,127,255}));
  connect(Sen_powlaw_1dat.port_b, bouB.ports[3]) annotation (Line(points={{12,48},
          {24,48},{24,-2.5},{60,-2.5}}, color={0,127,255}));
  connect(powlaw_2dat.port_b, Sen_powlaw_2dat.port_a)
    annotation (Line(points={{-20,10},{-8,10}},
                                              color={0,127,255}));
  connect(Sen_powlaw_2dat.port_b, bouB.ports[4]) annotation (Line(points={{12,10},
          {24,10},{24,-3.5},{60,-3.5}}, color={0,127,255}));
  connect(powlaw_M.port_b, Sen_powlaw_M.port_a)
    annotation (Line(points={{-20,-18},{-8,-18}},
                                                color={0,127,255}));
  connect(Sen_powlaw_M.port_b, bouB.ports[5]) annotation (Line(points={{12,-18},
          {24,-18},{24,-4.5},{60,-4.5}}, color={0,127,255}));
  connect(powlaw_V.port_b, Sen_powlaw_V.port_a)
    annotation (Line(points={{-20,-46},{-8,-46}},
                                                color={0,127,255}));
  connect(Sen_powlaw_V.port_b, bouB.ports[6]) annotation (Line(points={{12,-46},
          {28,-46},{28,-5.5},{60,-5.5}},
                                color={0,127,255}));
  connect(TabDat_M.port_b, Sen_tabdat_M.port_a)
    annotation (Line(points={{-20,-80},{-8,-80}},
                                                color={0,127,255}));
  connect(Sen_tabdat_M.port_b, bouB.ports[7]) annotation (Line(points={{12,-80},
          {28,-80},{28,-6.5},{60,-6.5}}, color={0,127,255}));
  connect(TabDat_V.port_b, Sen_tabdat_V.port_a)
    annotation (Line(points={{-20,-118},{-8,-118}}, color={0,127,255}));
  connect(Sen_tabdat_V.port_b, bouB.ports[8]) annotation (Line(points={{12,-118},
          {24,-118},{24,-7.5},{60,-7.5}}, color={0,127,255}));

  connect(Sen_ela.m_flow, m_flow_data[1]) annotation (Line(points={{2,117},{36,
          117},{36,83.25},{126,83.25}},
                                 color={0,0,127}));
  connect(Sen_ori.m_flow, m_flow_data[2]) annotation (Line(points={{2,89},{36,
          89},{36,85.75},{126,85.75}},
                                color={0,0,127}));
  connect(Sen_powlaw_M.m_flow, m_flow_data[3]) annotation (Line(points={{2,-7},{
          36,-7},{36,88.25},{126,88.25}},
                                    color={0,0,127}));
  connect(Sen_powlaw_V.m_flow, m_flow_data[4]) annotation (Line(points={{2,-35},
          {36,-35},{36,90.75},{126,90.75}},
                                     color={0,0,127}));
  connect(Sen_tabdat_M.m_flow, m_flow_data[5]) annotation (Line(points={{2,-69},
          {36,-69},{36,93.25},{126,93.25}},
                                     color={0,0,127}));
  connect(Sen_tabdat_V.m_flow, m_flow_data[6]) annotation (Line(points={{2,-107},
          {36,-107},{36,95.75},{126,95.75}},
                                      color={0,0,127}));
  connect(Sen_powlaw_2dat.m_flow, m_flow_data[7]) annotation (Line(points={{2,21},{
          36,21},{36,98.25},{126,98.25}},
                                        color={0,0,127}));
  connect(Sen_powlaw_1dat.m_flow, m_flow_data[8]) annotation (Line(points={{2,59},{
          36,59},{36,100.75},{126,100.75}},
                                        color={0,0,127}));
  connect(IntTestData.y, m_flow_testdata) annotation (Line(points={{-69,76},{-62,
          76},{-62,30},{36,30},{36,68},{126,68}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-140},
            {120,120}})),                                        Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-140,-140},{120,120}})),
    experiment(
      StopTime=500,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end OneWayFlow;
