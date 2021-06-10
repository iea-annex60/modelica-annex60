within IBPSA.Airflow.Multizone.BaseClasses;
partial model PartialPowerLawResistance_mflow
  "Partial model for mass flow resistance that uses the power law"
  extends IBPSA.Airflow.Multizone.BaseClasses.PartialOneWayFlowElement( final m_flow_nominal=rho_default*k*dp_turbulent);


  parameter Real m(min=0.5, max=1)
    "Flow exponent, m=0.5 for turbulent, m=1 for laminar";
  parameter Real k "Flow coefficient, k = m_flow/ dp^m";

protected
  constant Real gamma(min=1) = 1.5
    "Normalized flow rate where dphi(0)/dpi intersects phi(1)";

  parameter Real a = gamma
    "Polynomial coefficient for regularized implementation of flow resistance";
  parameter Real b = 1/8*m^2 - 3*gamma - 3/2*m + 35.0/8
    "Polynomial coefficient for regularized implementation of flow resistance";
  parameter Real c = -1/4*m^2 + 3*gamma + 5/2*m - 21.0/4
    "Polynomial coefficient for regularized implementation of flow resistance";
  parameter Real d = 1/8*m^2 - gamma - m + 15.0/8
    "Polynomial coefficient for regularized implementation of flow resistance";
equation

  m_flow = IBPSA.Airflow.Multizone.BaseClasses.powerLawFixedM(
    k=k,
    dp=dp,
    m=m,
    a=a,
    b=b,
    c=c,
    d=d,
    dp_turbulent=dp_turbulent);

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
end PartialPowerLawResistance_mflow;
