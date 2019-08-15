within IBPSA.Media.Refrigerants.R410A;
package R410A_IIR_P1_48_T233_340_Record "Refrigerant model for R410A using a hybrid approach with records"
  constant Modelica.Media.Interfaces.PartialTwoPhaseMedium.FluidConstants[1]
    refrigerantConstants(
     each chemicalFormula = "50% CH2F2 + 50% Pentafluorethan",
     each structureFormula = "50% Difluormethan + 50% CHF2CF3",
     each casRegistryNumber = "75-10-5 + 354-33-6",
     each iupacName = "R-410A,",
     each molarMass = 0.072585414240660,
     each criticalTemperature = 3.444943810810253e+02,
     each criticalPressure = 4.901264589893823e+06,
     each criticalMolarVolume = 1/6324,
     each normalBoilingPoint = 221.71,
     each triplePointTemperature = 200,
     each meltingPoint = 118.15,
     each acentricFactor = 0.296,
     each triplePointPressure = 29160,
     each dipoleMoment = 0,
     each hasCriticalData=true) "Thermodynamic constants for R410A";

  extends
    IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord(
  mediumName="R410A",
  substanceNames={"R410A"},
  singleState=false,
  SpecificEnthalpy(
    start=2.0e5,
    nominal=2.50e5,
    min=143.4e3,
    max=526.1e3),
  Density(
    start=500,
    nominal=750,
    min=5.1,
    max=1325),
  AbsolutePressure(
    start=2e5,
    nominal=5e5,
    min=1e5,
    max=48e5),
  Temperature(
    start=273.15,
    nominal=273.15,
    min=233.15,
    max=340),
  smoothModel=true,
  onePhase=false,
  ThermoStates=Choices.IndependentVariables.phX,
  fluidConstants=refrigerantConstants);

  redeclare record EoS
    "Record that contains fitting coefficients of the Helmholtz EoS"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R410A.EoS_IIR_P1_48_T233_340;
  annotation (Documentation(info="<html>
<p>This record contains fitting coefficients of the Helmholtz Equation of State.</p>
</html>"));
  end EoS;

  redeclare record BDSP
    "Record that contains fitting coefficients of the state properties at
    bubble and dew lines"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R410A.BDSP_IIR_P1_48_T233_340;
  annotation (Documentation(info="<html>
<p>This record contains fitting coefficients of the state properties at bubble and dew lines.</p>
</html>"));
  end BDSP;

  redeclare record TSP
    "Record that contains fitting coefficients of the state properties
    calculated with two independent state properties"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R410A.TSP_IIR_P1_48_T233_340;
  annotation (Documentation(info="<html>
<p>This record  contains fitting coefficients of the state properties calculated with two independent state properties.</p>
</html>"));
  end TSP;

  redeclare record SmoothTransition
    "Record that contains ranges to calculate a smooth transition between
    different regions"
      SpecificEnthalpy T_ph = 2.5;
      SpecificEntropy T_ps = 2.5;
      AbsolutePressure d_pT = 2.5;
      SpecificEnthalpy d_ph = 2.5;
      Real d_ps(unit="J/(Pa.K.kg)") =  50/(48e5-1e5);
      Real h_ps(unit="J/(Pa.K.kg)") = 50/(48e5-1e5);
  annotation (Documentation(info="<html>
<p>This record contains ranges to calculate a smooth transition between different regions.</p>
</html>"));
  end SmoothTransition;

  redeclare function extends dynamicViscosity
    "Calculates dynamic viscosity of refrigerant"

  /*The functional form of the dynamic viscosity is implented as presented in
    Geller et al. (2000), Viscosity of Mixed Refrigerants, R404A, R407C, R410A,
    and R507C. Eighth International Refrigeration Conference.
    Afterwards, the coefficients are adapted to the results obtained by the
    ExternalMedia libaray (i.e. CoolProp)
  */

protected
    SaturationProperties sat = setSat_T(state.T) "Saturation properties";
    Real phase_dT "Phase calculated by density and temperature";

    ThermodynamicState bubbleState "Thermodynamic state at bubble line";
    ThermodynamicState dewState "Thermodynamic state at dew line";
    Real quality "Vapour quality";

    Real etaZd "Dynamic viscosity for the limit of zero density";
    Real etaHd "Dynamic viscosity for the limit of high density";
    Real etaHdL
    "Dynamic viscosity for the limit of high density at bubble line";
    Real etaHdG "Dynamic viscosity for the limit of high density at dew line";
    Real etaL "Dynamic viscosity at dew line";
    Real etaG "Dynamic viscosity at bubble line";

  algorithm
    // Check phase
    if state.phase == 0 then
      phase_dT :=if not ((state.d < bubbleDensity(sat) and state.d > dewDensity(
        sat)) and state.T < fluidConstants[1].criticalTemperature) then 1 else 2;
    else
      phase_dT :=state.phase;
    end if;

    if (state.phase == 1 or phase_dT == 1) then
      // Calculate the dynamic viscosity near the limit of zero density
      etaZd := -2.695 + 5.850e-2*state.T - 2.129e-5*state.T^2;

      // Calculate the dynamic viscosity for limits of higher densities
      etaHd := 9.047e-3 + 5.784e-5*state.d^2 + 1.309e-7*state.d^3 -
        2.422e-10*state.d^4 + 9.424e-14*state.d^5 + 3.933e-17*state.d^6;

      // Calculate the final dynamic viscosity
      eta := (1.003684953*etaZd + 1.055260736*etaHd)*1e-6;
    else
      // Calculate properties at bubble and dew line
      bubbleState := setBubbleState(setSat_T(state.T));
      dewState := setDewState(setSat_T(state.T));
      quality := (bubbleState.d/state.d - 1)/(bubbleState.d/dewState.d - 1);

      // Calculate the dynamic viscosity near the limit of zero density
      etaZd := -2.695 + 5.850e-2*state.T - 2.129e-5*state.T^2;

      // Calculate the dynamic viscosity for limits of higher densities
      etaHdL := 9.047e-3 + 5.784e-5*bubbleState.d^2 + 1.309e-7*bubbleState.d^3 -
        2.422e-10*bubbleState.d^4 + 9.424e-14*bubbleState.d^5 +
        3.933e-17*bubbleState.d^6;
      etaHdG := 9.047e-3 + 5.784e-5*dewState.d^2 + 1.309e-7*dewState.d^3 -
        2.422e-10*dewState.d^4 + 9.424e-14*dewState.d^5 +
        3.933e-17*dewState.d^6;

      // Calculate the dynamic viscosity at bubble and dew line
      etaL := (1.003684953*etaZd + 1.055260736*etaHdL)*1e-6;
      etaG := (1.003684953*etaZd + 1.055260736*etaHdG)*1e-6;

      // Calculate the final dynamic viscosity
      eta := (quality/etaG + (1 - quality)/etaL)^(-1);
    end if;

  annotation (Documentation(info="<html>
<p>This function calculates dynamic viscosity of refrigerant.</p>
</html>"));
  end dynamicViscosity;

  redeclare function extends thermalConductivity
    "Calculates thermal conductivity of refrigerant"

  /*The functional form of the thermal conductivity is implented as presented
    in Geller et al. (2001), Thermal Conductivity of the Refrigerant Mixtures
    R404A, R407C, R410A, and R507A. International Journal of Thermophysics,
    Vol. 22, No. 4. Afterwards, the coefficients are adapted to the results
    obtained by the ExternalMedia libaray (i.e. CoolProp)
  */
protected
    SaturationProperties sat = setSat_T(state.T) "Saturation properties";
    Real phase_dT "Phase calculated by density and temperature";

    ThermodynamicState bubbleState "Thermodynamic state at bubble line";
    ThermodynamicState dewState "Thermodynamic state at dew line";
    Real quality "Vapour quality";

    Real lambdaIdg "Thermal conductivity for the limit of zero density";
    Real lambdaRes "Thermal conductivity for residual part";
    Real lambdaResL "Thermal conductivity for residual part at bubble line";
    Real lambdaResG "Thermal conductivity for residual part at dew line";
    Real lambdaL "Thermal conductivity at dew line";
    Real lambdaG "Thermal conductivity at bubble line";

  algorithm
    // Check phase
    if state.phase == 0 then
      phase_dT :=if not ((state.d < bubbleDensity(sat) and state.d > dewDensity(
        sat)) and state.T < fluidConstants[1].criticalTemperature) then 1 else 2;
    else
      phase_dT :=state.phase;
    end if;

    if (state.phase == 1 or phase_dT == 1) then
      // Calculate the thermal conductivity for the limit of zero density
      lambdaIdg := -8.872 + 7.41e-2*state.T;

      // Calculate the thermal conductivity for the residual part
      lambdaRes := 3.576e-2*state.d - 9.045e-6*state.d^2 + 4.343e-8*state.d^3
        - 3.705e-12*state.d^4;

      // Calculate the final thermal conductivity
      lambda := (lambdaIdg + 0.9994549202*lambdaRes)*1e-3;
    else
      // Calculate properties at bubble and dew line
      bubbleState := setBubbleState(setSat_T(state.T));
      dewState := setDewState(setSat_T(state.T));
      quality := (bubbleState.d/state.d - 1)/(bubbleState.d/dewState.d - 1);

      // Calculate the thermal conductivity for the limit of zero density
      lambdaIdg := -8.872 + 7.41e-2*state.T;

      // Calculate the thermal conductivity for the residual part
      lambdaResL := 3.576e-2*bubbleState.d - 9.045e-6*bubbleState.d^2 +
        4.343e-8*bubbleState.d^3 - 3.705e-12*bubbleState.d^4;
      lambdaResG := 3.576e-2*dewState.d - 9.045e-6*dewState.d^2 +
        4.343e-8*dewState.d^3 - 3.705e-12*dewState.d^4;

      // Calculate the thermal conductivity at bubble and dew line
      lambdaL := (lambdaIdg + 0.9994549202*lambdaResL)*1e-3;
      lambdaG := (lambdaIdg + 0.9994549202*lambdaResG)*1e-3;

      // Calculate the final thermal conductivity
      lambda := (quality/lambdaG + (1 - quality)/lambdaL)^(-1);
    end if;

  annotation (Documentation(info="<html>
<p>This function calculates thermal conductivity of refrigerant.</p>
</html>"));
  end thermalConductivity;

  redeclare function extends surfaceTension
    "Surface tension in two phase region of refrigerant"

  /*The functional form of the surface tension is implented as presented in
    Fröba and Leipertz (2003), Thermophysical Properties of the Refrigerant
    Mixtures R410A and R407C from Dynamic Light Scattering (DLS).
    International Journal ofThermophysics, Vol. 24, No. 5.
  */
protected
    Real tau = sat.Tsat/343.16 "Dimensionless temperature";

  algorithm
    sigma := (67.468*(1-tau)^1.26 * (1 - 0.051*(1-tau)^0.5 -
      0.193*(1-tau)))*1e-3;
  annotation (Documentation(info="<html>
<p>This function calculates the surface tension in two phase region of refrigerant.</p>
</html>"));
  end surfaceTension;

  annotation (Documentation(revisions="<html>
<ul>
<li>June 20, 2017, by Mirko Engelpracht:<br>First implementation (see <a href=\"https://github.com/RWTH-EBC/AixLib/issues/408\">issue 408</a>). </li>
<li>July 16, 2019, by Christian Vering</li>
</ul>
</html>", info="<html>
<p>
This package provides a refrigerant model for R410A using a hybrid approach
developed by Sangi et al.. The hybrid approach is implemented in
<a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMedium\">
IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMedium
</a>
and the refrigerant model is implemented by complete the template
<a href=\"modelica://IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord\">
IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord
</a>.
The fitting coefficients required in the template are saved in the package
<a href=\"modelica://IBPSA.DataBase.Media.Refrigerants.R410A\">
IBPSA.DataBase.Media.Refrigerants.R410A
</a>.
</p>
<h4>Assumptions and limitations</h4>
<p>
The implemented coefficients are fitted to external data by Engelpracht and
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
  <td><p>1</p></td>
  <td><p>48</p></td>
</tr>
<tr>
  <td><p>Temperature (T) in K</p></td>
  <td><p>233.15</p></td>
  <td><p>340.15</p></td>
</tr>
</table>
<p>
The reference point is defined as 200 kJ/kg and 1 kJ/kg/K, respectively, for
enthalpy and entropy for the saturated liquid at 273.15 K.
</p>
<h4>Assumptions and limitations</h4>
<p>
R410A is calculated as pseudo-pure fluid and, hence, only roughly valid
withing the two-phase region.
</p>
<h4>Validation</h4>
<p>
The model is validated by comparing results obtained from the example model
<a href=\"modelica://IBPSA.Media.Refrigerants.Examples.RefrigerantProperties\">
IBPSA.Media.Refrigerants.Examples.RefrigerantProperties
</a>
to external data (e.g. obtained from measurements or external media libraries).
</p>
<h4>References</h4>
<p>
Lemmon, E. W. (2003): Pseudo-Pure Fluid Equations of State for the Refrigerant
Blends R-410A, R-404A, R-507A, and R-407C. In: <i>International Journal of
Thermophysics 24 (4)</i>, S. 991–1006. DOI: 10.1023/A:1025048800563.
</p>
<p>
Geller, V. Z.; Bivens, D.; Yokozeki, A. (2000): Viscosity of Mixed
Refrigerants, R404A, R407C, R410A, and R507C. In: <i>International
refrigeration and air conditioning conference</i>. USA, S. 399–406. Online
available at http://docs.lib.purdue.edu/iracc/508.
</p>
<p>
Nabizadeh, H.; Mayinger, F. (1999): Viscosity of Gaseous R404A, R407C, R410A,
and R507. In: <i>International Journal of Thermophysics 20 (3)</i>, S.
777–790. DOI: 10.1007/978-1-4615-4777-8_1.
</p>
<p>
Geller, V. Z.; Nemzer, B. V.; Cheremnykh, U. V. (2001): Thermal Conductivity
of the Refrigerant Mixtures R404A, R407C, R410A, and R507A. In:
<i>International Journal of Thermophysics 22 (4)</i>, 1035–1043. DOI:
10.1023/A:1010691504352.
</p>
<p>
Fröba, A. P.; Leipertz, A. (2003): Thermophysical Properties of the
Refrigerant Mixtures R410A and R407C from Dynamic Light Scattering (DLS).
In: <i>International Journal ofThermophysics 24 (5)</i>, S. 1185–1206.
DOI: 10.1023/A:1026152331710.
</p>
<p>
Engelpracht, Mirko (2017): Development of modular and scalable simulation
models for heat pumps and chillers considering various refrigerants.
<i>Master Thesis</i>
</p>
</html>"));
end R410A_IIR_P1_48_T233_340_Record;
