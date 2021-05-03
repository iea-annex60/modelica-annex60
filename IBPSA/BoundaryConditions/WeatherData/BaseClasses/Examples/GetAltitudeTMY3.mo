within IBPSA.BoundaryConditions.WeatherData.BaseClasses.Examples;
model GetAltitudeTMY3 "Test model to get Altitude of TMY3"
  extends Modelica.Icons.Example;

  parameter String filNam = Modelica.Utilities.Files.loadResource(
  "modelica://IBPSA/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")
    "Name of weather data file";

  Modelica.SIunits.Length Altitude
    "Altitude of TMY3 location";



equation
  Altitude = IBPSA.BoundaryConditions.WeatherData.BaseClasses.getAltitudeLocationTMY3(
  filNam);


  annotation (
    Documentation(info="<html>
<p>
This example tests getting the time span of a TMY3 weather data file.
</p>
</html>",
revisions="<html>
<ul>
<li>
April 16, 2019, by Michael Wetter:<br/>
Removed call to get the absolute path of the file, corrected the <code>.mos</code>
file name and updated the documentation
</li>
<li>
April 15, 2019, by Ana Constantin:<br/>
First implementation.
</li>
</ul>
</html>"),
experiment(Tolerance=1e-6, StopTime=1.0),
__Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/BoundaryConditions/WeatherData/BaseClasses/Examples/GetTimeSpanTMY3.mos"
        "Simulate and plot"));
end GetAltitudeTMY3;
