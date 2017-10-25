within IBPSA.Fluid.FixedResistances.Validation.PlugFlowPipes;
model PlugFlowULg "Validation against data from Université de Liège"
  extends Modelica.Icons.Example;
  // R=((1/(2*pipe.kIns)*log((0.0603/2+pipe.dIns)/(0.0603/2)))+1/(5*(0.0603+2*pipe.dIns)))/Modelica.Constants.pi
  package Medium = IBPSA.Media.Water;

  parameter Modelica.SIunits.MassFlowRate m_flow_nominal=1
    "Nominal mass flow rate, used for regularization near zero flow";
  parameter Modelica.SIunits.Temperature T_start_in=pipeDataULg.T_start_in + 273.15
    "Initial temperature at pipe inlet";
  parameter Modelica.SIunits.Temperature T_start_out=pipeDataULg.T_start_out + 273.15
    "Initial temperature at pipe outlet";
  parameter Modelica.SIunits.SpecificHeatCapacity cp_default=
      Medium.specificHeatCapacityCp(state=sta_default)
    "Heat capacity of medium";
  parameter Medium.ThermodynamicState sta_default=Medium.setState_pTX(
      T=Medium.T_default,
      p=Medium.p_default,
      X=Medium.X_default) "Default medium state";

  Fluid.Sources.MassFlowSource_T WaterCityNetwork(
    redeclare package Medium = Medium,
    m_flow=1.245,
    use_m_flow_in=true,
    nPorts=1) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={70,0})));
  Fluid.HeatExchangers.Heater_T Boiler(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=0) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={44,0})));
  Fluid.Sources.Boundary_pT Sewer1(redeclare package Medium = Medium, nPorts=1)
    annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-110,0})));
  Fluid.Sensors.TemperatureTwoPort senTem_out(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0,
    T_start=T_start_out)
    annotation (Placement(transformation(extent={{-74,-10},{-94,10}})));
  Fluid.Sensors.TemperatureTwoPort senTem_in(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0,
    T_start=T_start_in)
    annotation (Placement(transformation(extent={{30,-10},{10,10}})));
  Modelica.Blocks.Sources.CombiTimeTable DataReader(table=pipeDataULg.data,
      extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
    annotation (Placement(transformation(extent={{0,-60},{20,-40}})));

  Modelica.Blocks.Math.UnitConversions.From_degC Tout
    "Ambient temperature in degrees"
    annotation (Placement(transformation(extent={{40,-88},{60,-68}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=295.15)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-30,70})));
  Modelica.Blocks.Math.UnitConversions.From_degC Tin
    "Input temperature into pipe"
    annotation (Placement(transformation(extent={{40,-60},{60,-40}})));

  replaceable Data.PipeDataULg151202 pipeDataULg constrainedby
    Data.BaseClasses.PipeDataULg "Measurement dataset from ULg"
    annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
  Modelica.Blocks.Math.Gain gain(k=1)
    annotation (Placement(transformation(extent={{52,-30},{72,-10}})));
  PlugFlowPipe pipe(
    redeclare package Medium = Medium,
    dh=0.05248,
    length=39,
    dIns(displayUnit="mm") = 0.013,
    kIns=0.04,
    m_flow_nominal=m_flow_nominal,
    thickness=3.9e-3,
    T_start_out=T_start_out,
    T_start_in=T_start_in,
    R=((1/(2*pipe.kIns)*log((0.0603/2 + pipe.dIns)/(0.0603/2))) + 1/(5*(0.0603
         + 2*pipe.dIns)))/Modelica.Constants.pi,
    nPorts=1,
    initDelay=true,
    m_flow_start=pipeDataULg.m_flowIni,
    cPip=500,
    rhoPip=8000)
    annotation (Placement(transformation(extent={{-20,-10},{-40,10}})));
  Fluid.Sensors.EnthalpyFlowRate senEntOut(redeclare package Medium = Medium,
      m_flow_nominal=m_flow_nominal)
    annotation (Placement(transformation(extent={{-46,-10},{-66,10}})));
  Modelica.Blocks.Math.Gain gain2(k=-1)
    annotation (Placement(transformation(extent={{20,40},{40,60}})));
  Modelica.Blocks.Math.MultiSum heatLossSim(nu=2)
    annotation (Placement(transformation(extent={{60,54},{72,66}})));
  Modelica.Blocks.Continuous.Integrator eneLosInt
    annotation (Placement(transformation(extent={{140,50},{160,70}})));
  Fluid.Sensors.EnthalpyFlowRate senEntIn(redeclare package Medium = Medium,
      m_flow_nominal=m_flow_nominal)
    annotation (Placement(transformation(extent={{4,-10},{-16,10}})));
  Modelica.Blocks.Math.Gain gain1(k=-1)
    annotation (Placement(transformation(extent={{74,-88},{94,-68}})));
  Modelica.Blocks.Math.MultiSum deltaT(nu=2)
    annotation (Placement(transformation(extent={{120,-80},{140,-60}})));
  Modelica.Blocks.Math.MultiProduct heatLossMeas(nu=2)
    annotation (Placement(transformation(extent={{160,-60},{180,-40}})));
  Modelica.Blocks.Math.Gain gain3(k=cp_default)
    annotation (Placement(transformation(extent={{120,-30},{140,-10}})));
  Modelica.Blocks.Math.Feedback heaLosDiff
    annotation (Placement(transformation(extent={{86,50},{106,70}})));
  Modelica.Blocks.Sources.Constant Tamb(k=273 + 18)
    "Ambient temperature in Kelvin";
equation
  connect(DataReader.y[3], Tout.u) annotation (Line(
      points={{21,-50},{32,-50},{32,-78},{38,-78}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(DataReader.y[5], Tin.u)
    annotation (Line(points={{21,-50},{29.5,-50},{38,-50}}, color={0,0,127}));
  connect(DataReader.y[1], gain.u) annotation (Line(points={{21,-50},{32,-50},{32,
          -20},{50,-20}}, color={0,0,127}));
  connect(senTem_in.port_a, Boiler.port_b)
    annotation (Line(points={{30,0},{30,0},{34,0}}, color={0,127,255}));
  connect(Boiler.port_a, WaterCityNetwork.ports[1])
    annotation (Line(points={{54,0},{54,0},{60,0}}, color={0,127,255}));
  connect(gain.y, WaterCityNetwork.m_flow_in) annotation (Line(points={{73,-20},
          {90,-20},{90,8},{80,8}}, color={0,0,127}));
  connect(Tin.y, Boiler.TSet) annotation (Line(points={{61,-50},{104,-50},{104,28},
          {62,28},{62,8},{56,8}},          color={0,0,127}));
  connect(Sewer1.ports[1], senTem_out.port_b)
    annotation (Line(points={{-100,0},{-94,0}}, color={0,127,255}));
  connect(senEntOut.H_flow, gain2.u) annotation (Line(points={{-56,11},{-56,26},
          {6,26},{6,50},{18,50}}, color={0,0,127}));
  connect(gain2.y, heatLossSim.u[1]) annotation (Line(points={{41,50},{50,50},{50,
          62.1},{60,62.1}}, color={0,0,127}));
  connect(senTem_out.port_a, senEntOut.port_b)
    annotation (Line(points={{-74,0},{-70,0},{-66,0}}, color={0,127,255}));
  connect(senEntOut.port_a, pipe.ports_b[1])
    annotation (Line(points={{-46,0},{-43,0},{-40,0}}, color={0,127,255}));
  connect(pipe.port_a, senEntIn.port_b)
    annotation (Line(points={{-20,0},{-16,0}}, color={0,127,255}));
  connect(senTem_in.port_b, senEntIn.port_a)
    annotation (Line(points={{10,0},{4,0}}, color={0,127,255}));
  connect(senEntIn.H_flow, heatLossSim.u[2]) annotation (Line(points={{-6,11},{-6,
          20},{52,20},{52,57.9},{60,57.9}}, color={0,0,127}));
  connect(fixedTemperature.port, pipe.heatPort)
    annotation (Line(points={{-30,60},{-30,10}}, color={191,0,0}));
  connect(Tout.y, gain1.u)
    annotation (Line(points={{61,-78},{72,-78}}, color={0,0,127}));
  connect(Tin.y, deltaT.u[1]) annotation (Line(points={{61,-50},{104,-50},{104,-66.5},
          {120,-66.5}},             color={0,0,127}));
  connect(gain1.y, deltaT.u[2]) annotation (Line(points={{95,-78},{104,-78},{104,
          -73.5},{120,-73.5}}, color={0,0,127}));
  connect(deltaT.y, heatLossMeas.u[1]) annotation (Line(points={{141.7,-70},{150,
          -70},{150,-46},{156,-46},{156,-46.5},{160,-46.5}},
                                                   color={0,0,127}));
  connect(gain.y, gain3.u)
    annotation (Line(points={{73,-20},{118,-20}},color={0,0,127}));
  connect(gain3.y, heatLossMeas.u[2]) annotation (Line(points={{141,-20},{148,-20},
          {148,-53.5},{160,-53.5}}, color={0,0,127}));
  connect(heatLossMeas.y, heaLosDiff.u2) annotation (Line(points={{181.7,-50},{200,
          -50},{200,40},{96,40},{96,52}},   color={0,0,127}));
  connect(heatLossSim.y, heaLosDiff.u1)
    annotation (Line(points={{73.02,60},{88,60}}, color={0,0,127}));
  connect(heaLosDiff.y, eneLosInt.u)
    annotation (Line(points={{105,60},{138,60}},          color={0,0,127}));
  annotation (
    Documentation(info="<html>
<p>
The example contains
experimental data from a real district heating network. 
</p>
<p>
This model compares the results with the original Modelica Standard Library pipes.
</p>
<p>The pipes' temperatures are not initialized. Therfore, results of
outflow temperature before approximately the first 10000 seconds should not be considered.
</p>
<h4>Test bench schematic</h4>
<p><img alt=\"Schematic of test rig at ULg\"
src=\"modelica://IBPSA/Resources/Images/Experimental/ULgTestBench.png\"/> </p>
<h4>Calibration</h4>
<p>
There are some uncertainties about the heat loss coefficient between pipe and surrounding air
as well as regarding the heat conductivity of the insulation material. 
With the <a href=\"modelica://IBPSA.Experimental.Pipe.Data.PipeDataULg150801\">given data</a>,
the length specific thermal resistance is <code>R = 2.164 </code>(mK/W), calculated as follows:
</p>
<p align=\"center\"style=\"font-style:italic;\">
R=((1/(2*pipe.kIns)*log((0.0603+2*pipe.dIns)/(0.0603)))+1/(5*(0.0603+2*pipe.dIns)))/Modelica.Constants.pi</p>
<p align=\"center\"style=\"font-style:italic;\">
U = 1/R = 0.462 W/(m K)</p>
</html>", revisions="<html>
<ul>
<li>November 24, 2016 by Bram van der Heijde:<br/>Add pipe thickness for wall capacity calculation and expand documentation section.</li>
<li>April 2, 2016 by Bram van der Heijde:<br/>Change thermal conductivity and put boundary condition in K.</li>
<li>Januar 26, 2016 by Carles Ribas:<br/>First implementation. </li>
</ul>
</html>"),
    experiment(StopTime=875, Tolerance=1e-006),
    __Dymola_Commands(file="Resources/Scripts/Dymola/Fluid/FixedResistances/Validation/PlugFlowPipes/PlugFlowULg.mos"
        "Simulate and plot"),
    Diagram(coordinateSystem(extent={{-140,-100},{220,100}})));
end PlugFlowULg;