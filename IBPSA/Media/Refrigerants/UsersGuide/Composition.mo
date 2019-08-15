within IBPSA.Media.Refrigerants.UsersGuide;
class Composition "Composition of the regrigerant library"
  extends Modelica.Icons.Information;

  annotation (Documentation(info="<html>
<p>The refrigerants&apos; library consists mainly of five packages. </p>
<ol>
<li><a href=\"IBPSA.Media.Refrigerants.Interfaces\">Interfaces</a><a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces\">:</a> Contains both templates to create new refrigerant models and partial models of the modeling approaches implemented in the library.</li>
<li><a href=\"IBPSA.Media.Refrigerants.DataBase\">DataBases</a><a href=\"modelica://IBPSA.DataBase.Media.Refrigerants\">:</a> Contains records with fitting coefficients used for, for example, different modeling approaches implemented in the library.</li>
<li><a href=\"IBPSA.Media.Refrigerants\">Refrigerants</a><a href=\"modelica://IBPSA.Media.Refrigerants.R134a\">:</a> Packages of different refrigerants which contain refrigerant models ready to use.</li>
<li><a href=\"IBPSA.Media.Refrigerants.Examples\">Examples</a><a href=\"modelica://IBPSA.Media.Refrigerants.Examples\">:</a> Contains example models to show the functionality of the refrigerant models.</li>
<li><a href=\"IBPSA.Media.Refrigerants.Examples\">Validation</a><a href=\"modelica://IBPSA.Media.Refrigerants.Validation\">:</a> Contains validation models to validate the modeling approaches implemented in the library.</li>
</ol>
<p>The ready to use models are provided in the following packages: </p>
<ul>
<li><a href=\"IBPSA.Media.Refrigerants.R134a\">R134a</a></li>
<li><a href=\"IBPSA.Media.Refrigerants.R290\">R290</a></li>
<li><a href=\"IBPSA.Media.Refrigerants.R410A\">R410A</a></li>
<li><a href=\"IBPSA.Media.Refrigerants.R32\">R32</a></li>
<li><a href=\"IBPSA.Media.Refrigerants.R744\">R744</a></li>
</ul>
</html>", revisions="<html>
<ul>
<li>October 14, 2017, by Mirko Engelpracht, Christian Vering:<br>First implementation (see <a href=\"https://github.com/RWTH-EBC/Aixlib/issues/408\">issue 408</a>). </li>
<li>July 7, 2019, by Bijan Sadjjadi</li>
<li>July 10, 2019, by Christian Vering</li>
</ul>
</html>"));
end Composition;
