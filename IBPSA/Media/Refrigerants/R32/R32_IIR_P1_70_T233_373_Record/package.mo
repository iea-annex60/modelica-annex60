within IBPSA.Media.Refrigerants.R32;
package R32_IIR_P1_70_T233_373_Record "Refrigerant model for R32 using a hybrid approach with records"

  /*Provide basic definitions of the refrigerant. Therefore, fill constants
    or parameters and may add new constants or parameters if needed. Moreover,
    provide references within the information of the package.
    */
    constant Modelica.Media.Interfaces.PartialTwoPhaseMedium.FluidConstants[1]
    refrigerantConstants(
    each chemicalFormula="CH2F2",
    each structureFormula="CH2F2",
    each casRegistryNumber="75-10-5",
    each iupacName="Difluoromethane",
    each molarMass=0.052024,
    each criticalTemperature=351.255,
    each criticalPressure=5782000,
    each criticalMolarVolume=0.052024/424,
    each triplePointTemperature=136.34,
    each triplePointPressure=47.9998938761,
    each normalBoilingPoint= 221,
    each meltingPoint= 137,
    each acentricFactor=0.2769,
    each dipoleMoment=0,
    each hasCriticalData=true) "Thermodynamic constants for R32";

    /*Provide basic information about the refrigerant. These basic information
    are the refrigerant name as well as the valid refrigerant limits in terms
    of specific enthalpy, density, absolute pressure and temperature.
  */
  extends
    IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMediumRecord(
  mediumName="R32",
  substanceNames={"R32"},
  singleState=false,
  SpecificEnthalpy(
    start=200000,
    nominal=200000,
    min=135338.010258576,
    max=627424.654013774),
  Density(
    start=500,
    nominal=350,
    min=1.68595438516567,
    max=1191.29136885037),
  AbsolutePressure(
    start=10e5,
    nominal=10e5,
    min=1e5,
    max=57.82e5),
  Temperature(
    start=268.37,
    nominal=333.15,
    min=233.15,
    max=373.15),
  smoothModel=true,
  onePhase=false,
  ThermoStates=Choices.IndependentVariables.phX,
  fluidConstants=refrigerantConstants);

    redeclare record EoS
    "Record that contains fitting coefficients of the Helmholtz EoS"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R32.EoS_IIR_P1_70_T233_373;
    annotation (Documentation(info=
                                 "<html>
<p>Record&nbsp;that&nbsp;contains&nbsp;fitting&nbsp;coefficients&nbsp;of&nbsp;the&nbsp;Helmholtz&nbsp;EoS.</p>
</html>"));
    end EoS;

  redeclare record BDSP
    "Record that contains fitting coefficients of the state properties at
    bubble and dew lines"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R32.BDSP_IIR_P1_70_T233_373;
  annotation (Documentation(info="<html>
<p>Record&nbsp;that&nbsp;contains&nbsp;fitting&nbsp;coefficients&nbsp;of&nbsp;the&nbsp;state&nbsp;properties&nbsp;at bubble&nbsp;and&nbsp;dew&nbsp;lines.</p>
</html>"));
  end BDSP;

  redeclare record TSP
    "Record that contains fitting coefficients of the state properties
    calculated with two independent state properties"
    extends
      IBPSA.Media.Refrigerants.DataBase.Refrigerants.R32.TSP_IIR_P1_70_T233_373;
  annotation (Documentation(info="<html>
<p>Record&nbsp;that&nbsp;contains&nbsp;fitting&nbsp;coefficients&nbsp;of&nbsp;the&nbsp;state&nbsp;properties calculated&nbsp;with&nbsp;two&nbsp;independent&nbsp;state&nbsp;properties.</p>
</html>"));
  end TSP;

  redeclare record SmoothTransition "Record that contains ranges to calculate a 
  smooth transition between different regions"
    SpecificEnthalpy T_ph=2.5;
    SpecificEntropy T_ps=2.5;
    AbsolutePressure d_pT=2.5;
    SpecificEnthalpy d_ph=2.5;
    Real d_ps(unit="J/(Pa.K.kg)") = 25/(39.5e5 - 1e5);
    Real h_ps(unit="J/(Pa.K.kg)") = 25/(39.5e5 - 1e5);
  annotation (Documentation(info="<html>
<p>Record&nbsp;that&nbsp;contains&nbsp;ranges&nbsp;to&nbsp;calculate&nbsp;a&nbsp;smooth&nbsp;transition&nbsp;between different&nbsp;regions.</p>
</html>"));
  end SmoothTransition;

redeclare function extends dynamicViscosity
    "Calculates dynamic viscosity of refrigerant"

    /*The functional form of the dynamic viscosity is implented as presented in
 Huber et al. (1999), Transport propertiers of refrigerants R32, R125, R134a and R125 + R32 mixtures in and beyond the critical region */

protected
    SaturationProperties sat = setSat_T(state.T) "Saturation properties";
    Real phase_dT "Phase calculated by density and temperature";

    ThermodynamicState bubbleState "Thermodynamic state at bubble line";
    ThermodynamicState dewState "Thermodynamic state at dew line";
    Real quality "Vapour quality";
    Real Tred = state.T/262.6 "Reduced temperature for lower density terms";
    Real omegaEta "Reduced effective collision cross section";
    Real dilgasEta "dilute-gas viscosity";
    Real M = 52 "molar Mass in g/mol";
    Real resEta "Residual portion of the viscosity";
    Real resEtaL "Residual portion of the viscosity";
    Real resEtaG "Residual portion of the viscosity";
    Real etaL;
    Real etaG;
    Real dcrit = 424 / 0.0520 * 1e-3 "critical density";
    Real dmol = state.d / 0.0520 * 1e-3 "molar density in mol/L";
    Real dmol_bub "molar density at bubble line in mol/L";
    Real dmol_dew "molar density at dew line in mol/L";

algorithm
    // Check phase
    if state.phase == 0 then
      phase_dT :=if not ((state.d < bubbleDensity(sat) and state.d > dewDensity(sat)) and state.T < fluidConstants[1].criticalTemperature) then 1 else 2;
    else
      phase_dT :=state.phase;
    end if;
    if (state.phase == 1 or phase_dT == 1) then
    // Calculate the dynamic visocity near the limit of zero density
      if abs(Tred)<1E-20 then
        Tred := 1E-20;
      end if;
    omegaEta :=1.16145/Tred^(0.14874) + 0.52487*exp((-0.77320)*Tred) + 2.16178*exp(-2.43787*Tred);
    dilgasEta :=26.69579e-9*sqrt(Tred*M)/(omegaEta*0.4282^2);
    resEta :=10^(-7)*exp(0.3745666e1 + (-0.2593564e4)/state.T)*(exp((0.3861782e1 + 0.1807046e5/state.T^(1.5))
        *dmol^(0.1) + (dmol/dcrit - 1)*sqrt(dmol)*(0.1478233e1 + (-0.5001262e3)/state.T + 0.5455701e5/state.T
        ^2)) - 1);
    eta :=dilgasEta + resEta;

    else
  // Calculate properties at bubble and dew line
      bubbleState := setBubbleState(setSat_T(state.T));
      dewState := setDewState(setSat_T(state.T));
      quality := (bubbleState.d/state.d - 1)/(bubbleState.d/dewState.d - 1);

      // Calculate the dynamic visocity near the limit of zero density
      if abs(Tred)<1E-20 then
        Tred := 1E-20;
      end if;
      dmol_bub :=bubbleState.d/0.0520*1e-3;
      dmol_dew :=dewState.d/0.0520*1e-3;
      omegaEta :=1.16145/Tred^(0.14874) + 0.52487*exp((-0.77320)*Tred) + 2.16178*exp(-2.43787*Tred);
      dilgasEta :=26.69579e-9*sqrt(Tred*M)/(omegaEta*0.4282^2);

      resEtaL :=10^(-7)*exp(0.3745666e1 + (-0.2593564e4)/state.T)*(exp((0.3861782e1 + 0.1807046e5/state.T^(1.5))
        *dmol_bub^(0.1) + (dmol_bub/dcrit - 1)*sqrt(dmol_bub)*(0.1478233e1 + (-0.5001262e3)/state.T + 0.5455701e5
        /state.T^2)) - 1);
      resEtaG :=10^(-7)*exp(0.3745666e1 + (-0.2593564e4)/state.T)*(exp((0.3861782e1 + 0.1807046e5/state.T^(1.5))
        *dmol_dew^(0.1) + (dmol_dew/dcrit - 1)*sqrt(dmol_dew)*(0.1478233e1 + (-0.5001262e3)/state.T + 0.5455701e5
        /state.T^2)) - 1);

      etaL :=dilgasEta + resEtaL;
      etaG :=dilgasEta + resEtaG;

      eta := (quality/etaG + (1 - quality)/etaL)^(-1);
    end if;

  annotation (Documentation(info="<html>
<p>Calculates&nbsp;dynamic&nbsp;viscosity&nbsp;of&nbsp;refrigerant.</p>
</html>"));
end dynamicViscosity;

  redeclare function extends thermalConductivity
  "Calculates thermal conductivity of refrigerant"
    /*The functional form of the thermal conductivity is implented as presented in
    Huber et al. (1999), Transport propertiers of refrigerants R32, R125, R134a and R125 + R32 mixtures in and beyond the critical region*/

protected
      SaturationProperties sat = setSat_T(state.T) "Saturation properties";
      Integer phase_dT "Phase calculated by density and temperature";

      ThermodynamicState bubbleState "Thermodynamic state at bubble line";
      ThermodynamicState dewState "Thermodynamic state at dew line";
      Real quality "Vapour quality";
      Real lambdaIdg "Thermal conductivity for the limit of zero density";
      Real lambdaRes "Thermal conductivity for residual part";
      Real lambdaResL "Thermal conductivity for residual part at bubble line";
      Real lambdaResG "Thermal conductivity for residual part at dew line";
      Real lambdaG;
      Real lambdaL;
      Real dmol = state.d / 0.0520 * 1e-3 "molar density in mol/L";
      Real dmol_bub "molar density at bubble line in mol/L";
      Real dmol_dew "molar density at dew line in mol/L";

  algorithm
      // Check phase
      if state.phase == 0 then
          phase_dT :=if not ((state.d < bubbleDensity(sat) and state.d > dewDensity(sat)) and state.T < fluidConstants[1].criticalTemperature) then 1 else 2;
        else
          phase_dT :=state.phase;
        end if;
      if (state.phase == 1 or phase_dT == 1) then
        lambdaIdg :=0.124363e-1 + (-0.628224e-3)*state.T + 0.199947e-6*state.T^2;
        lambdaRes :=(0.793100e-2 + (-0.129997e-4)*state.T)*dmol + ((-0.373196e-3) + 0.189222e-5*state.T)*
        dmol^2 + (0.195280e-4 + (-0.116382e-6)*state.T)*dmol^3 + ((-0.267308e-6) + 0.299505e-8*state.T)*dmol^
        4;
        lambda :=lambdaIdg + lambdaRes;

      else
        bubbleState := setBubbleState(setSat_T(state.T));
        dewState := setDewState(setSat_T(state.T));
        quality := (bubbleState.d/state.d - 1)/(bubbleState.d/dewState.d - 1);

        // Calculate the dynamic visocity near the limit of zero density
        dmol_bub :=bubbleState.d/0.0520*1e-3;
        dmol_dew :=dewState.d/0.0520*1e-3;

        lambdaIdg :=0.124363e-1 + (-0.628224e-3)*state.T + 0.199947e-6*state.T^2;
        lambdaResL :=(0.793100e-2 + (-0.129997e-4)*state.T)*dmol_bub + ((-0.373196e-3) + 0.189222e-5*state.T)
        *dmol_bub^2 + (0.195280e-4 + (-0.116382e-6)*state.T)*dmol_bub^3 + ((-0.267308e-6) + 0.299505e-8*
        state.T)*dmol_bub^4;
        lambdaResG :=(0.793100e-2 + (-0.129997e-4)*state.T)*dmol_dew + ((-0.373196e-3) + 0.189222e-5*state.T)
        *dmol_dew^2 + (0.195280e-4 + (-0.116382e-6)*state.T)*dmol_dew^3 + ((-0.267308e-6) + 0.299505e-8*
        state.T)*dmol_dew^4;

        lambdaL :=lambdaIdg + lambdaResL;
        lambdaG :=lambdaIdg + lambdaResG;

        // Calculate the final dynamic visocity
        lambda := (quality/lambdaG + (1 - quality)/lambdaL)^(-1);
      end if;

  annotation (Documentation(info="<html>
<p>Calculates&nbsp;thermal&nbsp;conductivity&nbsp;of&nbsp;refrigerant.</p>
</html>"));
  end thermalConductivity;

  redeclare function extends surfaceTension
  "Surface tension in two phase region of refrigerant"

    /*The functional form of the surface tension is implemented as presented in
    Fröba and Leipertz (2003), Thermophysical Properties of the Refrigerant
    Mixtures R410A and R407C from Dynamic Light Scattering (DLS).
    International Journal ofThermophysics, Vol. 24, No. 5.
  */
protected
  Real x = abs(1- sat.Tsat/343.16) "Dimensionless temperature";

  algorithm
    sigma :=0.07147*x^(1.246);

  annotation (Documentation(info="<html>
<p>Surface&nbsp;tension&nbsp;in&nbsp;two&nbsp;phase&nbsp;region&nbsp;of&nbsp;refrigerant.</p>
</html>"));
  end surfaceTension;

annotation (Documentation(info="<html>
<p>This package provides a refrigerant model for R32 using a hybrid approach developed by Sangi et al.. The hybrid approach is implemented in IBPSA.Media.Refrigerants.Interfaces.PartialHybridTwoPhaseMedium and the refrigerant model is implemented by complete the template IBPSA.Media.Refrigerants.Interfaces.TemplateHybridTwoPhaseMediumRecord . The fitting coefficients required in the template are saved in the package IBPSA.DataBase.Media.Refrigerants.R32 . </p>
</html>", revisions="<html>
<ul>
<li>December 10, 2018, by Stephan Goebel: First implementation (see <a href=\"https://github.com/RWTH-EBC/Aixlib/issues/665\">issue 665</a>)</li>
<li>July 16, 2019, by Christian Vering</li>
</ul>
</html>
"));
end R32_IIR_P1_70_T233_373_Record;
