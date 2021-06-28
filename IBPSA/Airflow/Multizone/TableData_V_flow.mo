within IBPSA.Airflow.Multizone;
model TableData_V_flow
  "Volume flow(y-axis) vs Pressure(x-axis) cubic spline fit model based on table data, with last two points linearly interpolated"

  extends IBPSA.Airflow.Multizone.BaseClasses.PartialOneWayFlowElement(
    V_flow = IBPSA.Airflow.Multizone.BaseClasses.flowElementData(u=dp,xd=xd,yd=yd,d=d),
    final m_flow_nominal=min(abs(table[:,2]))/rho_default);

parameter Real table[:,:] "Table of volume flow rate in m3/s (second column) as a function of pressure difference in Pa (first column)";

protected
  parameter   Real[:] xd=table[:,1] "X-axis support points";
  parameter   Real[size(xd, 1)] yd=table[:,2] "Y-axis support points";
  parameter   Real[size(xd, 1)] d(each fixed=false)
                                                   "Derivatives at the support points";

initial equation
  d =IBPSA.Utilities.Math.Functions.splineDerivatives(
    x=xd,
    y=yd,
    ensureMonotonicity=true);

  annotation (
    defaultComponentName="tabdat_V",
    Documentation(info="<html>
<p>
This model describes the one-directional pressure driven air flow through an 
opening based on fixed tabular input describing the relation between mass flow rate
and the pressure difference over the component.
</p>
<p>
<img src=\"modelica://IBPSA/Resources/Images/equations/equation-tGaXKOnB.png\" alt=\"V_flow = f(dp)\"/>
</p>
<p>
<i>dp = the pressure difference over the flow element</i>
</p>
<p>
<i>V_flow = the volume flow through the element (positive from A-&gt;B)</i>
</p>
<p>
Based on the table input, a cubic hermite spline is constructed between all point 
except for the second to last and last point. These point are connected linearly.
</p>
<p>
The constructed curve is the direct relation between dp and V_flow.
</p>
<p>
A similar model is also used in the CONTAM software (Dols and Walton, 2015).
</p>
<h4>References</h4>
<ul>
<li>
<b>W. S. Dols and B. J. Polidoro</b>,<b>2015</b>. <i>CONTAM User Guide and Program Documentation Version 3.2</i>, 
National Institute of Standards and Technology, NIST TN 1887, Sep. 2015. doi: 
<a href=\"https://doi.org/10.6028/NIST.TN.1887\">10.6028/NIST.TN.1887</a>. 
</li>
</ul>
</html>", revisions="<html>
<ul>
<li>Apr 6, 2021, 2020, by Klaas De Jonge (UGent):<br/>
First implementation</li>
</ul>
</html>
"), Icon(graphics={
        Rectangle(
          extent={{-48,80},{52,-80}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-98,6},{-48,-6}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,127,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{52,6},{102,-6}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,127,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-22,78},{2,-78}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{2,78},{28,-78}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{0,13},{0,-13}},
          lineColor={0,0,127},
          origin={15,0},
          rotation=90,
          textString="V_flow"),
        Rectangle(
          extent={{-22,78},{28,58}},
          lineColor={28,108,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Forward),
        Line(points={{2,78},{2,-78}}, color={28,108,200})}));
end TableData_V_flow;
