within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.BorefieldData;
record ExampleBorefieldData
  "Example definition of a borefield data record"
  extends
    IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.BorefieldData.Template(
      filDat=IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.FillingData.Bentonite(),
      soiDat=IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.SoilData.SandStone(),
      conDat=IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.ConfigurationData.ExampleConfigurationData());
  annotation (
defaultComponentPrefixes="parameter",
defaultComponentName="borFieDat",
Documentation(
info="<html>
<p>This record presents an example on how to define BorefieldData records
using the template in
<a href=\"modelica://IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.BorefieldData.Template\">
IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.BorefieldData.Template</a>.</p>
</html>",
revisions="<html>
<ul>
<li>
July 15, 2018, by Michael Wetter:<br/>
Revised implementation, added <code>defaultComponentPrefixes</code> and
<code>defaultComponentName</code>.
</li>
<li>
June 28, 2018, by Damien Picard:<br/>
First implementation.
</li>
</ul>
</html>"));
end ExampleBorefieldData;
