within IBPSA.BoundaryConditions.WeatherData.BaseClasses;
function getAltitudeLocationTMY3 "Gets the altitude from TMY3 file"
  extends Modelica.Icons.Function;
 input String filNam "Name of weather data file"
 annotation (Dialog(
        loadSelector(filter="Weather files (*.mos)", caption=
            "Select weather file")));
 output Modelica.SIunits.Length Altitude "Altitude of TMY3 location";
protected
 Integer nexInd "Next index, used for error handling";
 String element "String representation of the returned element";
algorithm
  element :=
    IBPSA.BoundaryConditions.WeatherData.BaseClasses.getLastHeaderElementTMY3(
      filNam=filNam,
      start="#LOCATION",
      name = "Altitude");
   (nexInd, Altitude) :=Modelica.Utilities.Strings.Advanced.scanReal(
    string=element,
    startIndex=1,
    unsigned=false);

   // Check if altitude is valid
   assert(Altitude > 0,
       "Wrong value for Altitude. Received Altitude = " +
       String(Altitude));

  annotation (Documentation(info="<html>
This function returns the time zone of the TMY3 weather data file.
</html>", revisions="<html>
<ul>
<li>
May 2, 2021, by Ettore Zanetti:<br/>
First implementation.
</li>
</ul>
</html>"));
end getAltitudeLocationTMY3;
