within Modelica_Fluid;
package Types "Common types for fluid models"

  annotation (preferedView="info",
              Documentation(info="<html>
 
</html>"));

  type HydraulicConductance =Real (
      final quantity="HydraulicConductance",
      final unit="kg/(s.Pa)");
  type HydraulicResistance =Real (
      final quantity="HydraulicResistance",
      final unit="Pa.s/kg");

  type FrictionTypes = enumeration(
      ConstantLaminar "constant laminar flow",
      ConstantTurbulent "constant turbulent",
      DetailedFriction "detailed friction model")
    "Enumeration to define the pressure loss equations due to friction";

  type CrossSectionTypes = enumeration(
      Circular "circular",
      Rectangular "rectangular",
      General "general")
    "Enumeration to define the geometric cross section of pipes";

  type Dynamics = enumeration(
      DynamicFreeInitial
        "DynamicFreeInitial -- Dynamic balance, Initial guess value",
      FixedInitial "FixedInitial -- Dynamic balance, Initial value fixed",
      SteadyStateInitial
        "SteadyStateInitial -- Dynamic balance, Steady state initial with guess value", 

      SteadyState "SteadyState -- Steady state balance, Initial guess value")
    "Enumeration to define definition of balance equations" 
  annotation (Documentation(info="<html>
<p>
Enumeratioin to define the formulation of balance equations
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.Dynamics.</b></th><th><b>Meaning</b></th></tr>
<tr><td>DynamicFreeInitial</td><td>Dynamic balance, Initial guess value</td></tr>
 
<tr><td>FixedInitial</td><td>Dynamic balance, Initial value fixed</td></tr>
 
<tr><td>SteadyStateInitial</td><td>Dynamic balance, Steady state initial with guess value</td></tr>
 
<tr><td>SteadyState</td><td>Steady state balance, Initial guess value</td></tr>
</table>
</html>"));

  type CvTypes = enumeration(
      Av "Av (metric) flow coefficient",
      Kv "Kv (metric) flow coefficient",
      Cv "Cv (US) flow coefficient",
      OpPoint "Av defined by operating point")
    "Enumeration to define the choice of valve flow coefficient";

  type PortFlowDirection = enumeration(
      Entering "Fluid flow is only entering",
      Leaving "Fluid flow is only leaving",
      Bidirectional "No restrictions on fluid flow (flow reversal possible)")
    "Enumeration to define whether flow reversal is allowed";

  type ModelStructure = enumeration(
      av_vb "av_vb: port_a - volume - flow model - volume - port_b",
      a_v_b "a_v_b: port_a - flow model - volume - flow model - port_b",
      av_b "av_b: port_a - volume - flow model - port_b",
      a_vb "a_vb: port_a - flow model - volume - port_b")
    "Enumeration with choices for model structure in distributed pipe model";
end Types;
