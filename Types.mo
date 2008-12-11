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
      Dynamic
        "Dynamic -- Differential equations with accumulation of mass, energy and momentum as applicable", 

      SteadyStateMomentum
        "SteadyStateMomentum -- Steady state momentum equations, accumulation of mass and energy", 

      SteadyStateMass
        "SteadyStateMass -- No variation of mass, steady state momentum, accumulation of energy", 

      SteadyState
        "SteadyState -- Algebraic model, no accumulation of mass, energy or momentum")
    "Enumeration to define dynamics options" 
  annotation (Documentation(info="<html>
<p>
Enumeration type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.Dynamics.</b></th><th><b>Meaning</b></th></tr>
<tr><td>Dynamic</td>
    <td>Differential equations with accumulation of mass, energy and momentum as applicable</td></tr>
 
<tr><td>SteadyStateMomentum</td>
    <td>Steady state momentum equations, accumulation of mass and energy</td></tr>
 
<tr><td>SteadyStateHoldup</td>
    <td>No variation of mass, steady state momentum, accumulation of energy</td></tr>
 
<tr><td>SteadyState</td>
    <td>Algebraic model, no accumulation of mass, energy or momentum</td></tr>
</table>
</html>"));

  type Init = enumeration(
      InitialValues "InitialValues -- Initial values for p, T or h, X",
      NoInit "NoInit -- No initial conditions (guess values for p, T or h, X)", 

      SteadyStateMomentum "SteadyStateMomentum: Steady state momentum",
      SteadyStateHydraulic
        "SteadyStateHydraulic -- Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X", 

      SteadyState "SteadyState -- Steady state (guess values for p, T or h, X)")
    "Enumeration to define initialization options" 
  annotation (Documentation(info="<html>
<p>
Integer type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.Init.</b></th><th><b>Meaning</b></th></tr>
<tr><td>InitialValues</td>
    <td>Initial values for p, T or h, X</td></tr>
 
<tr><td>NoInit</td>
    <td>No initial conditions (guess values for p, T or h, X)</td></tr>
 
<tr><td>SteadyStateMomentum</td>
    <td>Steady state momentum</td></tr>
 
<tr><td>SteadyStateHydraulic</td>
    <td>Hydraulic steady state (der(p)=0), guess value for p, 
        initial values for T or h, X</td></tr>
 
<tr><td>SteadyState</td>
    <td>Steady state (guess values for p, T or h, X)</td></tr>
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
