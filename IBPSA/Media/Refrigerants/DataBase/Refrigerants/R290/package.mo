within IBPSA.Media.Refrigerants.DataBase.Refrigerants;
package R290 "Package provides records for R290"
  extends Modelica.Icons.VariantsPackage;



  annotation (Documentation(info="<html>
<p>
This package provides records with fitting coefficients for the refrigerant
R290. The records are inherited from the base data definitions provided in
IBPSA.DataBase.Media.Refrigerants and the fitting coefficients are used for
the refrigerant model provided in IBPSA.Media.Refrigerants.
</p>
<p>
For detailed information of the <b>base data definitions</b>, please
checkout the following records:
</p>
<ul>
<li><a href=\"modelica://IBPSA.DataBase.Media.Refrigerants.HelmholtzEquationOfStateBaseDateDefinition\">
IBPSA.DataBase.Media.Refrigerants.HelmholtzEquationOfStateBaseDateDefinition
</a></li>
<li><a href=\"modelica://IBPSA.DataBase.Media.Refrigerants.BubbleDewStatePropertiesBaseDataDefinition\">
IBPSA.DataBase.Media.Refrigerants.BubbleDewStatePropertiesBaseDataDefinition
</a></li>
<li><a href=\"modelica://IBPSA.DataBase.Media.Refrigerants.ThermodynamicStatePropertiesBaseDataDefinition\">
IBPSA.DataBase.Media.Refrigerants.ThermodynamicStatePropertiesBaseDataDefinition
</a></li>
</ul>
<p>
For detailed information of the <b>refrigerant models using the fitting
coefficients</b>, please checkout the following packages:
</p>
<ul>
<li><a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord\">
IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord
</a></li>
<li><a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord\">
IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord
</a></li>
</ul>
</html>", revisions="<html>
<ul>
<li>June 10, 2017, by Mirko Engelpracht:<br>First implementation (see <a href=\"https://github.com/RWTH-EBC/AixLib/issues/408\">issue 408</a>). </li>
<li>July 15, 2019, by Christian Vering</li>
</ul>
</html>"));
end R290;
