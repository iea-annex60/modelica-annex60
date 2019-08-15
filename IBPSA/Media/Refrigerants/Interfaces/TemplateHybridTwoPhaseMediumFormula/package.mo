within IBPSA.Media.Refrigerants.Interfaces;
partial package TemplateHybridTwoPhaseMediumFormula "Template for media models using a hybrid approach without records"
  constant Modelica.Media.Interfaces.PartialTwoPhaseMedium.FluidConstants[1]
    refrigerantConstants(
      each chemicalFormula = "CXHY",
      each structureFormula = "CXHY",
      each casRegistryNumber = "xx-xx-x",
      each iupacName = "name",
      each molarMass = 1,
      each criticalTemperature = 1,
      each criticalPressure = 1,
      each criticalMolarVolume = 1,
      each normalBoilingPoint = 1,
      each triplePointTemperature = 1,
      each meltingPoint = 1,
      each acentricFactor = 1,
      each triplePointPressure = 1,
      each dipoleMoment = 1,
      each hasCriticalData=true)
    "Thermodynamic constants for refrigerant";

  extends
    IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula(
    mediumName="Name",
    substanceNames={"Name"},
    singleState=false,
    SpecificEnthalpy(
      start=1.0e5,
      nominal=1.0e5,
      min=50e3,
      max=1000e3),
    Density(
      start=500,
      nominal=500,
      min=0.5,
      max=100),
    AbsolutePressure(
      start=1e5,
      nominal=5e5,
      min=0.5e5,
      max=30e5),
    Temperature(
      start=273.15,
      nominal=273.15,
      min=200.15,
      max=423.15),
    smoothModel=true,
    onePhase=false,
    ThermoStates=.Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.phX,
    fluidConstants=refrigerantConstants);

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

  redeclare function extends f_Idg
    "Dimensionless Helmholtz energy (Ideal gas contribution f_Idg)"
  algorithm

  end f_Idg;

  redeclare function extends f_Res
    "Dimensionless Helmholtz energy (Residual part alpha_r)"
  algorithm

  end f_Res;

  redeclare function extends t_fIdg_t
    "Short form for tau*(dalpha_0/dtau)_delta=const"
  algorithm

  end t_fIdg_t;

  redeclare function extends tt_fIdg_tt
    "Short form for tau*tau*(ddalpha_0/(dtau*dtau))_delta=const"
  algorithm

  end tt_fIdg_tt;

  redeclare function extends t_fRes_t
    "Short form for tau*(dalpha_r/dtau)_delta=const"
  algorithm

  end t_fRes_t;

  redeclare function extends tt_fRes_tt
    "Short form for tau*tau*(ddalpha_r/(dtau*dtau))_delta=const"
  algorithm

  end tt_fRes_tt;

  redeclare function extends d_fRes_d
    "Short form for delta*(dalpha_r/(ddelta))_tau=const"
  algorithm

  end d_fRes_d;

  redeclare function extends dd_fRes_dd
    "Short form for delta*delta(ddalpha_r/(ddelta*delta))_tau=const"
  algorithm

  end dd_fRes_dd;

  redeclare function extends td_fRes_td
    "Short form for tau*delta*(ddalpha_r/(dtau*ddelta))"
  algorithm

  end td_fRes_td;

  redeclare function extends ttt_fIdg_ttt
    "Short form for tau*tau*tau*(dddalpha_0/(dtau*dtau*dtau))_delta=const"
  algorithm

  end ttt_fIdg_ttt;

  redeclare function extends ttt_fRes_ttt
    "Short form for tau*tau*tau*(dddalpha_r/(dtau*dtau*dtau))_delta=const"
  algorithm

  end ttt_fRes_ttt;

  redeclare function extends ddd_fRes_ddd
    "Short form for delta*delta*delta*
    (dddalpha_r/(ddelta*ddelta*ddelta))_tau=const"
  algorithm

  end ddd_fRes_ddd;

  redeclare function extends tdd_fRes_tdd
    "Short form for tau*delta*delta*(dddalpha_r/(dtau*ddelta*ddelta))"
  algorithm

  end tdd_fRes_tdd;

  redeclare function extends ttd_fRes_ttd
    "Short form for tau*tau*delta*(dddalpha_r/(dtau*dtau*ddelta))"
  algorithm

  end ttd_fRes_ttd;

  redeclare function extends saturationPressure
    "Saturation pressure of refrigerant (Ancillary equation)"
protected
    Real OM = (1 - T/fluidConstants[1].criticalTemperature);

  algorithm
    if T>fluidConstants[1].criticalTemperature then
     p := fluidConstants[1].criticalPressure;
    elseif T<fluidConstants[1].triplePointTemperature then
     p := fluidConstants[1].triplePointPressure;
    else

    end if;

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end saturationPressure;

  redeclare function extends saturationTemperature
    "Saturation temperature of refrigerant (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end saturationTemperature;

  redeclare replaceable function extends bubbleDensity
    "Boiling curve specific density of refrigerant (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end bubbleDensity;

  redeclare function extends dewDensity
    "Dew curve specific density of refrigerant (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end dewDensity;

  redeclare function extends bubbleEnthalpy
    "Boiling curve specific enthalpy of refrigerant (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end bubbleEnthalpy;

  redeclare function extends dewEnthalpy
    "Dew curve specific enthalpy of refrigerant (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end dewEnthalpy;

  redeclare function extends bubbleEntropy
    "Boiling curve specific entropy of refrigerant (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end bubbleEntropy;

  redeclare function extends dewEntropy
    "Dew curve specific entropy of propane (Ancillary equation)"
protected
    Real x = 0;

  algorithm

    annotation(smoothOrder = 2,
               Inline=false,
               LateInline=true);
  end dewEntropy;

  redeclare function temperature_ph
    "Calculates temperature as function of pressure and specific enthalpy"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input FixedPhase phase = 0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output Temperature T "Temperature";

protected
    SmoothTransition st;
    SpecificEnthalpy dh = st.T_ph;
    SpecificEnthalpy h_dew = 0;
    SpecificEnthalpy h_bubble = 0;
    Real x1 = 0;
    Real y1 = 0;
    Real T1 = 0;
    Real x2 = 0;
    Real y2 = 0;
    Real T2 = 0;

  algorithm

  annotation(derivative(noDerivative=phase)=temperature_ph_der,
             inverse(h=specificEnthalpy_pT(p=p,T=T,phase=phase)));
  end temperature_ph;

  redeclare function temperature_ps
    "Calculates temperature as function of pressure and specific entroy"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input FixedPhase phase = 0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output Temperature T "Temperature";

protected
    SmoothTransition st;
    SpecificEntropy ds = st.T_ps;
    SpecificEntropy s_dew = 0;
    SpecificEntropy s_bubble = 0;
    Real x1 = 0;
    Real y1 = 0;
    Real T1 = 0;
    Real x2 = 0;
    Real y2 = 0;
    Real T2 = 0;

  algorithm

  annotation(derivative(noDerivative=phase)=temperature_ps_der);
  end temperature_ps;

  redeclare function density_pT
    "Computes density as a function of pressure and temperature"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input FixedPhase phase = 0
      "2 for two-phase, 1 for one-phase, 0 if not known";
    output Density d "Density";

protected
    SmoothTransition st;
    AbsolutePressure dp = st.d_pT;
    SaturationProperties sat = setSat_T(T=T);
    Real x1 = 0;
    Real y1 = 0;
    Real d1 = 0;
    Real x2 = 0;
    Real y2 = 0;
    Real d2 = 0;

  algorithm

  annotation(derivative(noDerivative=phase)=density_pT_der,
             inverse(p=pressure_dT(d=d,T=T,phase=phase)));
  end density_pT;

  redeclare function extends dynamicViscosity
    "Calculates dynamic viscosity of refrigerant"

  algorithm

  end dynamicViscosity;

  redeclare function extends thermalConductivity
    "Calculates thermal conductivity of refrigerant"

  algorithm

  end thermalConductivity;

  redeclare function extends surfaceTension
    "Surface tension in two phase region of refrigerant"

  algorithm

  end surfaceTension;

  annotation (Documentation(revisions="<html>
<ul>
<li>June 6, 2017, by Mirko Engelpracht:<br>First implementation (see <a href=\"https://github.com/RWTH-EBC/Aixlib/issues/408\">issue 408</a>). </li>
<li>July 16, 2019, by Christian Vering</li>
</ul>
</html>",
        info="<html>
<p>
This package is a <b>template</b> for <b>new refrigerant</b> models using a
hybrid approach based on the &quot;Fast_Propane&quot; model developed by
Sangi et al. (for detailed information, please checkout
<a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula\">
IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumFormula
</a>
). For a new refrigerant model just make a copy of this package, remove the
&quot;partial&quot; keyword from the package and provide the information that
is requested in the comments of the Modelica code. The following
<b>information is requested</b>:
</p>
<ol>
<li><u>Fluid constants:</u> Provide basic information of the refrigerant and
add the reference.</li>
<li><u>Basic information:</u> Provide basic information like the refrigerant
name and its valid fluid limits in terms of different thermodynamic state
properties.</li>
<li><u>Base properties:</u> Provide information about the refrigerant&apos;s
base properties like the relation between specific enthalpy, specific internal
energy, pressure and density.</li>
<li><u>Helmholtz EoS:</u> Provide the formulas of the Helmholtz equation
of state.</li>
<li><u>Saturation state properties:</u> Provide the formulas of the
saturation properties.</li>
<li><u>Further thermodynamic properties:</u> Provide formulas for further
thermodynamic properties like the thermal conductivity or surface tension.</li>
<li><u>Fitted state properties:</u> Provide the formulas of ther further
state properties.</li>
</ol>
<p>
A refrigerant package inherits from <b>PartialHybridTwoPhaseMediumFormula</b>
and provides the equations for the refrigerant. Moreover, the
PartialHybridTwoPhaseMedium package inherits from <b>PartialMedium</b> and,
therefore, the details of this package are described in
<a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide
</a>.
</p>
<h4>References</h4>
<p>
Sangi, Roozbeh; Jahangiri, Pooyan; Klasing, Freerk; Streblow, Rita;
M&uuml;ller, Dirk (2014): <a href=\"http://dx.doi.org/10.3384/ecp14096\">
A Medium Model for the Refrigerant Propane for Fast and Accurate Dynamic
Simulations</a>. In: <i>The 10th International Modelica Conference</i>. Lund,
Sweden, March 10-12, 2014: Link&ouml;ping University Electronic Press
(Link&ouml;ping Electronic Conference Proceedings), S. 1271&ndash;1275
</p>
</html>"));
end TemplateHybridTwoPhaseMediumFormula;
