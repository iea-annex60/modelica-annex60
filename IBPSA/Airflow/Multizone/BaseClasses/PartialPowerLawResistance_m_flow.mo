within IBPSA.Airflow.Multizone.BaseClasses;
partial model PartialPowerLawResistance_m_flow
  "Partial model for mass flow resistance that uses the power law for computing mass flow rate"
  extends IBPSA.Airflow.Multizone.BaseClasses.PartialOneWayFlowElement(
    m_flow = IBPSA.Airflow.Multizone.BaseClasses.powerLawFixedM(
      k=k,
      dp=dp,
      m=m,
      a=a,
      b=b,
      c=c,
      d=d,
      dp_turbulent=dp_turbulent),
    final m_flow_nominal=k*dp_turbulent,
    final m_flow_small=1E-4*abs(m_flow_nominal));
  extends IBPSA.Airflow.Multizone.BaseClasses.PowerLawResistanceParameters;


  parameter Real k "Flow coefficient, k = m_flow/ dp^m";

  annotation (
    Documentation(info="<html>
<p>This model describes the mass flow rate and pressure difference relation in the form </p>
<p align=\"center\" style=\"font-style:italic;\">
m_flow = k * dp^m,
</p>
<p>where <i>k</i> and <i>m</i> are parameters. For turbulent flow, set <i>m=1/2</i> and for laminar flow, set <i>m=1</i>. </p>
<p>The model is used as a base for the interzonal air flow models. </p>
</html>",
revisions="<html>
<ul>
<li>
Apr 6, 2021, by Klaas De Jonge (UGent):<br/>
Power law resistance describing mass flow. This model is the same as <a href=\"modelica://IBPSA.Airflow.Multizone.BaseClasses.PowerLawResistance\"> IBPSA.Airflow.Multizone.BaseClasses.PowerLawResistance</a> but with the mass flow rate as a function of dP in stead of volume flow rate
</li>
</ul>
</html>

"));
end PartialPowerLawResistance_m_flow;
