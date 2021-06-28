within IBPSA.Airflow.Multizone.BaseClasses.Examples;
model WindPressureProfile
  "Test model for wind pressure profile function"
  extends Modelica.Icons.Example;


  parameter Real table[:,:]=
    [0,0.4; 45,0.1; 90,-0.3; 135,-0.35; 180,-0.2; 225,-0.35; 270,-0.3; 315,0.1; 360,0.4]
    "Table data";
  Modelica.SIunits.Angle  alpha "Wind incidence angle (0: normal to wall)";
  Real Cp "Wind pressure coefficient";


protected
  Real Radtable[:,:] = [Modelica.Constants.D2R*table[:,1],table[:,2]]
  "Transforms table with input in degrees to radians";

  //Extend table to ensure correct derivatives at 0 and 360.
  Real prevPoint[1,2] = [Radtable[size(table, 1)-1, 1] - (2*Modelica.Constants.pi), Radtable[size(table, 1)-1, 2]]
  "Second to last point of the input table";
  Real nextPoint[1,2] = [Radtable[2, 1] + (2*Modelica.Constants.pi), Radtable[2, 2]]
  "Second point of the input table";
  Real exTable[:,:]=[prevPoint;Radtable;nextPoint]
  "Add additional points to the transformed table";

  //Arguments for windPressureProfile function
  Real[:] xd=exTable[:,1] "Support points x-value";
  Real[size(xd, 1)] yd=exTable[:,2] "Support points y-value";
  Real[size(xd, 1)] d=IBPSA.Utilities.Math.Functions.splineDerivatives(x=xd,y=yd,ensureMonotonicity=false) "Spline derivative values at the support points";

  Modelica.Blocks.Sources.Ramp ramp(
    duration=500,
    height=3*360,
    offset=-360)
    "Ramp model generating a singal from -360 to 720";


initial equation
  //continuity at 0 and 360 must be assured
  assert(table[1,2]<>0 or table[end,2]<>360, "First and last point in the table must be 0 and 360", level = AssertionLevel.error);

equation
   alpha=Modelica.Constants.D2R*ramp.y;
   Cp =IBPSA.Airflow.Multizone.BaseClasses.windPressureProfile(u=alpha, xd=xd,yd=yd,d=d);

  annotation (
experiment(
      StopTime=500,
      Tolerance=1e-06),
  __Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Airflow/Multizone/BaseClasses/Examples/WindPressureProfile.mos"
        "Simulate and plot"), Documentation(info="<html>
<p>
This examples demonstrates the
<a href=\"modelica://IBPSA.Airflow.Multizone.BaseClasses.windPressureProfile\">
IBPSA.Airflow.Multizone.BaseClasses.windPressureProfile</a>
function.
</p>
</html>", revisions="<html>
<ul>
<li>
Apr 6, 2021, by Klaas De Jonge (UGent):<br/>
First implementation
</li>
</ul>
</html>
"));
end WindPressureProfile;
