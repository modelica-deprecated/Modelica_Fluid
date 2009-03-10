within Modelica_Fluid;
package Media
  replaceable partial model BaseProperties
    "Computes the basic thermodynamic properties corresponding to a thermodynamic state"
    import IV = ModelicaNew.Media.Interfaces.Types.IndependentVariables;
    replaceable package Medium = ModelicaNew.Media.Interfaces.GenericMedium;
    constant Integer nXi = Medium.nS - 1 "Number of independent mass fractions";
    InputAbsolutePressure p "Absolute pressure of medium";
    InputSpecificEnthalpy h "Specific enthalpy of medium";
    InputMassFraction Xi[nXi](start=Medium.reference_X[1:nXi])
      "Independent mass fractions";
    Medium.Density d "Density of medium";
    Medium.Temperature T "Temperature of medium";
    Medium.MassFraction[Medium.nS] X(start=Medium.reference_X)
      "Mass fractions (= (component mass)/total mass  m_i/m)";
    Medium.SpecificInternalEnergy u "Specific internal energy of medium";
    // SpecificHeatCapacity R "Gas constant (of mixture if applicable)";
    // MolarMass MM "Molar mass (of mixture or single fluid)";
    Medium.ThermodynamicState state
      "thermodynamic state record for optional functions";
    parameter Boolean preferredMediumStates=false
      "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
      annotation (Evaluate=true, Dialog(tab="Advanced"));
    annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={0,0,255}), Text(
            extent={{-152,164},{152,102}},
            textString="%name",
            lineColor={0,0,255})}));

    // Local connector definition, used for equation balancing check
    connector InputAbsolutePressure = input SI.AbsolutePressure
      "Pressure as input signal connector";
    connector InputSpecificEnthalpy = input SI.SpecificEnthalpy
      "Specific enthalpy as input signal connector";
    connector InputMassFraction = input SI.MassFraction
      "Mass fraction as input signal connector";

  equation
    // Calculation of thermodynamic state based on ThermoState
    if Medium.ThermoStates == IV.pTX then
      state = Medium.setState_pTX(p,T,X);
      h = Medium.specificEnthalpy(state);
    elseif Medium.ThermoStates == IV.pTY then
      state = Medium.setState_phX(p,h,X);
      T = Medium.temperature(state);
    else
      assert(false, "Unsupported choice of thermo states");
    end if;

    // Calculation of properties
    d = Medium.density(state);
    u = Medium.specificEnthalpy(state)-p/d;
    // MM = Medium.molarMass(state);
    // R= Medium.R;
    // Total quantities
    X = cat(1,Xi,{1-sum(Xi)});

  /*
  if standardOrderComponents then
    Xi = X[1:nXi];
 
      if fixedX then
        X = reference_X;
      end if;
      if reducedX and not fixedX then
        X[nX] = 1 - sum(Xi);
      end if;
      for i in 1:nX loop
        assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X[" +
               String(i) + "] = " + String(X[i]) + "of substance "
               + substanceNames[i] + "\nof medium " + mediumName + " is not in the range 0..1");
      end for;
 
  end if;
*/
    assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" +
      Medium.mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
    annotation (Documentation(info="<html>
<p>
Model <b>BaseProperties</b> is a model within package <b>PartialMedium</b>
and contains the <b>declarations</b> of the minimum number of
variables that every medium model is supposed to support.
A specific medium inherits from model <b>BaseProperties</b> and provides
the equations for the basic properties. </p>
<p>
The BaseProperties model contains the following <b>7+nXi variables</b>
(nXi is the number of independent mass fractions defined in package
PartialMedium):
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>Variable</b></td>
      <td valign=\"top\"><b>Unit</b></td>
      <td valign=\"top\"><b>Description</b></td></tr>
  <tr><td valign=\"top\">T</td>
      <td valign=\"top\">K</td>
      <td valign=\"top\">temperature</td></tr>
  <tr><td valign=\"top\">p</td>
      <td valign=\"top\">Pa</td>
      <td valign=\"top\">absolute pressure</td></tr>
  <tr><td valign=\"top\">d</td>
      <td valign=\"top\">kg/m^3</td>
      <td valign=\"top\">density</td></tr>
  <tr><td valign=\"top\">h</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific enthalpy</td></tr>
  <tr><td valign=\"top\">u</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific internal energy</td></tr>
  <tr><td valign=\"top\">Xi[nXi]</td>
      <td valign=\"top\">kg/kg</td>
      <td valign=\"top\">independent mass fractions m_i/m</td></tr>
  <tr><td valign=\"top\">R</td>
      <td valign=\"top\">J/kg.K</td>
      <td valign=\"top\">gas constant</td></tr>
  <tr><td valign=\"top\">M</td>
      <td valign=\"top\">kg/mol</td>
      <td valign=\"top\">molar mass</td></tr>
</table>
<p>
In order to implement an actual medium model, one can extend from this
base model and add <b>5 equations</b> that provide relations among 
these variables. Equations will also have to be added in order to
set all the variables within the ThermodynamicState record state.</p>
<p>
If standardOrderComponents=true, the full composition vector X[nX] 
is determined by the equations contained in this base class, depending 
on the independent mass fraction vector Xi[nXi]. </p>
<p>Additional <b>2 + nXi</b> equations will have to be provided 
when using the BaseProperties model, in order to fully specify the 
thermodynamic conditions. The input connector qualifier applied to 
p, h, and nXi indirectly declares the number of missing equations,
permitting advanced equation balance checking by Modelica tools.
Please note that this doesn't mean that the additional equations 
should be connection equations, nor that exactly those variables
should be supplied, in order to complete the model.
For further information, see the Modelica.Media User's guide, and 
Section 4.7 (Balanced Models) of the Modelica 3.0 specification. </p>
</html>"));
  end BaseProperties;
end Media;
