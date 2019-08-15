within IBPSA.Media.Refrigerants.R290;
package R290_IIR_P05_30_T263_343_Record "Refrigerant model for R290 using a hybrid approach with records developed
     by Sangi et al."
  constant Modelica.Media.Interfaces.PartialTwoPhaseMedium.FluidConstants[1]
    refrigerantConstants(
     each chemicalFormula = "C3H8",
     each structureFormula = "C3H8",
     each casRegistryNumber = "74-98-6",
     each iupacName = "Propane",
     each molarMass = 0.04409562,
     each criticalTemperature = 369.89,
     each criticalPressure = 4.2512e6,
     each criticalMolarVolume = 1/(5e3),
     each normalBoilingPoint = 231.036,
     each triplePointTemperature = 85.525,
     each meltingPoint = 85.45,
     each acentricFactor = 0.153,
     each triplePointPressure = 0.00017,
     each dipoleMoment = 0.1,
     each hasCriticalData=true) "Thermodynamic constants for Propane";

  extends
    IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord(
  mediumName="Propane",
  substanceNames={"Propane"},
  singleState=false,
  SpecificEnthalpy(
    start=2.057e5,
    nominal=2.057e5,
    min=177e3,
    max=576e3),
  Density(
    start=300,
    nominal=529,
    min=0.77,
    max=547),
  AbsolutePressure(
    start=4.7446e5,
    nominal=5e5,
    min=0.5e5,
    max=30e5),
  Temperature(
    start=273.15,
    nominal=333.15,
    min=263.15,
    max=343.15),
  smoothModel=true,
  onePhase=false,
  ThermoStates=Choices.IndependentVariables.phX,
  fluidConstants=refrigerantConstants);

  redeclare record EoS
    "Record that contains fitting coefficients of the Helmholtz EoS"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R290.EoS_IIR_P05_30_T263_343;
  end EoS;

  redeclare record BDSP
    "Record that contains fitting coefficients of the state properties at
    bubble and dew lines"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R290.BDSP_IIR_P05_30_T263_343;
  end BDSP;

  redeclare record TSP
    "Record that contains fitting coefficients of the state properties
    calculated with two independent state properties"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R290.TSP_IIR_P05_30_T263_343;
  end TSP;

  redeclare record SmoothTransition
    "Record that contains ranges to calculate a smooth transition between
    different regions"
      SpecificEnthalpy T_ph = 2.5;
      SpecificEntropy T_ps = 2.5;
      AbsolutePressure d_pT = 2.5;
      SpecificEnthalpy d_ph = 2.5;
      Real d_ps(unit="J/(Pa.K.kg)") =  25/(30e5-0.5e5);
      Real h_ps(unit="J/(Pa.K.kg)") = 50/(30e5-0.5e5);
  end SmoothTransition;

  redeclare function extends dynamicViscosity
    "Calculates dynamic viscosity of refrigerant"

protected
    Real tv[:] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0,
                  2.0, 3.0, 4.0, 1.0, 2.0};
    Real dv[:] = {1.0, 2.0, 3.0, 13.0, 12.0, 16.0, 0.0, 18.0,
                  20.0, 13.0, 4.0, 0.0, 1.0};
    Real nv[:] = {-0.7548580e-1, 0.7607150, -0.1665680, 0.1627612e-5,
                  0.1443764e-4, -0.2759086e-6, -0.1032756, -0.2498159e-7,
                  0.4069891e-8, -0.1513434e-5, 0.2591327e-2, 0.5650076,
                  0.1207253};

    Real d_crit = MM/fluidConstants[1].criticalMolarVolume "Critical density";
    Real MM = fluidConstants[1].molarMass "Molar mass";

    ThermodynamicState dewState = setDewState(setSat_T(state.T))
    "Thermodynamic state at dew line";
    ThermodynamicState bubbleState = setBubbleState(setSat_T(state.T))
    "Thermodynamic state at bubble line";
    Real dr "Reduced density for higher density terms";
    Real drL "Reduced density for higher density terms at bubble line";
    Real drG "Reduced density for higher density terms at dew line";
    Real etaL
    "Dynamic viscosity for the limit of high density at bubble line";
    Real etaG "Dynamic viscosity for the limit of high density at dew line";
    Real Hc = 17.1045;
    Real Tr = state.T/fluidConstants[1].criticalTemperature
    "Reduced temperature for higher density terms";

    SaturationProperties sat = setSat_T(state.T);
    Real phase_dT = if not ((state.d < bubbleDensity(sat) and state.d >
      dewDensity(sat)) and state.T < fluidConstants[1].criticalTemperature)
      then 1 else 2;
    Real quality;

  algorithm
    if state.phase==1 or phase_dT==1 then
      eta := 0;
      dr := state.d/d_crit;
      for i in 1:11 loop
          eta := eta + nv[i]*Tr^tv[i]*dr^dv[i];
      end for;
      for i in 12:13 loop
          eta := eta + exp(-dr*dr/2)*nv[i]*Tr^tv[i]*dr^dv[i];
      end for;
      eta := (exp(eta) - 1)*Hc/1e6;
    elseif state.phase==2 or phase_dT==2 then
      etaL := 0;
      etaG := 0;
      drG := dewState.d/d_crit;
      drL := bubbleState.d/d_crit;
      quality :=if state.phase == 2 then (bubbleState.d/state.d - 1)/(
      bubbleState.d/dewState.d - 1) else 1;
      for i in 1:11 loop
          etaL := etaL + nv[i]*Tr^tv[i]*drL^dv[i];
          etaG := etaG + nv[i]*Tr^tv[i]*drG^dv[i];
      end for;
      for i in 12:13 loop
          etaL := etaL + exp(-drL*drL/2)*nv[i]*Tr^tv[i]*drL^dv[i];
          etaG := etaG + exp(-drG*drG/2)*nv[i]*Tr^tv[i]*drG^dv[i];
      end for;
      etaL := (exp(etaL) - 1)*Hc/1e6;
      etaG := (exp(etaG) - 1)*Hc/1e6;
      eta := (quality/etaG + (1 - quality)/etaL)^(-1);
    end if;
  end dynamicViscosity;

  redeclare function extends thermalConductivity
    "Calculates thermal conductivity of refrigerant"

protected
    Real B1[:] = {-3.51153e-2,1.70890e-1,-1.47688e-1,5.19283e-2,-6.18662e-3};
    Real B2[:] = {4.69195e-2,-1.48616e-1,1.32457e-1,-4.85636e-2,6.60414e-3};
    Real C[:] = {3.66486e-4,-2.21696e-3,2.64213e+0};
    Real A[:] = {-1.24778e-3,8.16371e-3,1.99374e-2};

    Real d_crit = MM/fluidConstants[1].criticalMolarVolume;
    Real MM = fluidConstants[1].molarMass;

    Real delta "Reduced density for higher density terms";
    Real deltaL "Reduced density for higher density terms at bubble line";
    Real deltaG "Reduced density for higher density terms at dew line";
    Real tau = fluidConstants[1].criticalTemperature/state.T
    "Reduced temperature for higher density terms";

    Real lambda0 = A[1]+A[2]/tau+A[3]/(tau^2);
    Real lambdar = 0;
    Real lambdarL = 0;
    Real lambdarG = 0;
    Real lambdaL = 0 "Thermal conductivity at bubble line";
    Real lambdaG = 0 "Thermal conductivity at dew line";

    SaturationProperties sat = setSat_T(state.T);
    Real phase_dT = if not ((state.d < bubbleDensity(sat) and state.d >
      dewDensity(sat)) and state.T < fluidConstants[1].criticalTemperature)
      then 1 else 2;
    Real quality;

  algorithm
    delta :=state.d/d_crit;
    if state.phase==1 or phase_dT==1 then
      for i in 1:5 loop
          lambdar := lambdar + (B1[i] + B2[i]/tau)*delta^i;
      end for;
      lambda := (lambda0 + lambdar + (C[1]/(C[2] + abs(1.0/tau - 1.0))*
        exp(-(C[3]*(delta - 1.0))^2)));
    elseif state.phase==2 or phase_dT==2 then
      deltaL :=bubbleDensity(setSat_T(state.T))/d_crit;
      deltaG :=dewDensity(setSat_T(state.T))/d_crit;
      quality :=if state.phase == 2 then (bubbleDensity(setSat_T(state.T))/
        state.d - 1)/(bubbleDensity(setSat_T(state.T))/dewDensity(setSat_T(
        state.T)) - 1) else 1;
      for i in 1:5 loop
          lambdarL := lambdarL + (B1[i] + B2[i]/tau)*deltaL^i;
          lambdarG := lambdarG + (B1[i] + B2[i]/tau)*deltaG^i;
      end for;
      lambdaL := (lambda0 + lambdarL + (C[1]/(C[2] + abs(1.0/tau - 1.0))*
        exp(-(C[3]*(deltaL - 1.0))^2)));
      lambdaG := (lambda0 + lambdarG + (C[1]/(C[2] + abs(1.0/tau - 1.0))*
        exp(-(C[3]*(deltaG - 1.0))^2)));
      lambda := (quality/lambdaG + (1 - quality)/lambdaL)^(-1);
    end if;
  end thermalConductivity;

  redeclare function extends surfaceTension
    "Surface tension in two phase region of refrigerant"

  algorithm
    sigma := 1e-3*55.817*(1-sat.Tsat/369.85)^1.266;
  end surfaceTension;

  annotation (Documentation(revisions="<html>
<ul>
<li>June 12, 2017, by Mirko Engelpracht:<br>First implementation (see <a href=\"https://github.com/RWTH-EBC/Aixlib/issues/408\">issue 408</a>). </li>
<li>July 16, 2019, by Christian Vering</li>
</ul>
</html>", info="<html>
<p>
This package provides a refrigerant model for R290 using a hybrid approach
developed by Sangi et al.. The hybrid approach is implemented in
<a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord\">
IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord
</a>
and the refrigerant model is implemented by complete the template
<a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord\">
IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord
</a>.
The fitting coefficients required in the template are saved in the package
<a href=\"modelica://IBPSA.DataBase.Media.Refrigerants.R290\">
IBPSA.DataBase.Media.Refrigerants.R290</a>.
</p>
<h4>Assumptions and limitations</h4>
<p>
The implemented coefficients are fitted to external data by Sangi et al. and
are valid within the following range:<br />
</p>
<table summary=\"Range of validiry\" cellspacing=\"0\" cellpadding=\"2\"
border=\"1\" width=\"30%\" style=\"border-collapse:collapse;\">
<tr>
  <td><p>Parameter</p></td>
  <td><p>Minimum Value</p></td>
  <td><p>Maximum Value</p></td>
</tr>
<tr>
  <td><p>Pressure (p) in bar</p></td>
  <td><p>0.5</p></td>
  <td><p>30</p></td>
</tr>
<tr>
  <td><p>Temperature (T) in K</p></td>
  <td><p>263.15</p></td>
  <td><p>343.15</p></td>
</tr>
</table>
<h4>Validation</h4>
<p>
Sangi et al. validated their model by comparing it to results obtained from
the Helmholtz equation of state. They found out that relative error of the
refrigerant model compared to HelmholtzMedia (Thorade and Saadat, 2012) is
close to zero.
</p>
<h4>References</h4>
<p>
Thorade, Matthis; Saadat, Ali (2012):
<a href=\"http://www.ep.liu.se/ecp/076/006/ecp12076006.pdf\">
HelmholtzMedia - A fluid properties library</a>. In: <i>Proceedings of the
9th International Modelica Conference</i>; September 3-5; 2012; Munich;
Germany. Link&ouml;ping University Electronic Press, S. 63&ndash;70.
</p>
<p>
Sangi, Roozbeh; Jahangiri, Pooyan; Klasing, Freerk; Streblow, Rita;
M&uuml;ller, Dirk (2014): <a href=\"http://dx.doi.org/10.3384/ecp14096\">
A Medium Model for the Refrigerant Propane for Fast and Accurate Dynamic
Simulations</a>. In: <i>The 10th International Modelica Conference</i>. Lund,
Sweden, March 10-12, 2014: Link&ouml;ping University Electronic Press
(Link&ouml;ping Electronic Conference Proceedings), S. 1271&ndash;1275
</p>
<p>
Klasing,Freerk: A New Design for Direct Exchange Geothermal Heat Pumps -
Modeling, Simulation and Exergy Analysis. <i>Master thesis</i>
</p>
<p>
Scalabrin, G.; Marchi, P.; Span, R. (2006): A Reference Multiparameter
Viscosity Equation for Propane with an Optimized Functional Form. In: <i>J.
Phys. Chem. Ref. Data, Vol. 35, No. 3</i>, S. 1415-1442
</p>
</html>"));
end R290_IIR_P05_30_T263_343_Record;
