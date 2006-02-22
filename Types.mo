package Types "Common types for fluid models" 
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
Package <b>Types</b> contains common type definitions of the Modelica_Fluid
library.
</p>
</html>"));
  
  type HydraulicConductance =Real (
      final quantity="HydraulicConductance",
      final unit="kg/(s.Pa)");
  type HydraulicResistance =Real (
      final quantity="HydraulicResistance",
      final unit="Pa.s/kg");
  
  package FrictionTypes 
    "Type, constants and menu choices to define the pressure loss equations due to friction, as temporary solution until enumerations are available" 
    
    extends Modelica.Icons.Enumeration;
    constant Integer ConstantLaminar=1;
    constant Integer ConstantTurbulent=2;
    constant Integer DetailedFriction=3;
    type Temp 
      "Temporary type of FrictionTypes with choices for menus (until enumerations are available)" 
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar 
            "ConstantLaminar \"dp = k*m_flow\"",
          choice=Modelica_Fluid.Types.FrictionTypes.ConstantTurbulent 
            "ConstantTurbulent \"dp = k*m_flow^2\"",
          choice=Modelica_Fluid.Types.FrictionTypes.DetailedFriction 
            "DetailedFriction \"dp = f(Re,delta,rho,L,D,nu)\""));
    end Temp;
  end FrictionTypes;
  
  package CrossSectionTypes 
    "Type, constants and menu choices to define the geometric cross section of pipes, as temporary solution until enumerations are available" 
    
    annotation (preferedView="text");
    
    extends Modelica.Icons.Enumeration;
    constant Integer Circular=1;
    constant Integer Rectangular=2;
    constant Integer General=3;
    type Temp 
      "Temporary type of CrossSectionTypes with choices for menus (until enumerations are available)" 
      
      extends Integer;
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.CrossSectionTypes.Circular 
            "Circular cross section",
          choice=Modelica_Fluid.Types.CrossSectionTypes.Rectangular 
            "Rectangular cross section",
          choice=Modelica_Fluid.Types.CrossSectionTypes.General 
            "General cross section"));
    end Temp;
  end CrossSectionTypes;
  
  package InitTypes "Obsolete (will be removed), use Types.Init" 
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit = 0 "No explicit initial conditions";
    constant Integer InitialValues = 1 
      "Initial conditions specified by initial values";
    constant Integer SteadyState = 2 "Full steady-state";
    constant Integer SteadyStateHydraulic = 3 
      "Hydraulic steady state, initial values of other state variables given";
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=0,max=3);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.Init.NoInit 
            "NoInit (No explicit initial conditions)",
          choice=Modelica_Fluid.Types.Init.InitialValues 
            "InitialValues (Initial conditions specified by initial values)",
          choice=Modelica_Fluid.Types.Init.SteadyState 
            "SteadyState (Full steady state)",
          choice=Modelica_Fluid.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (Hydraulic steady state, initial values of other states given)"));
    end Temp;
  end InitTypes;
  
  package Init 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    
    annotation (Documentation(info="<html>
 
</html>"));
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit = 1 
      "No initial conditions (guess values for p, T or h, X)";
    constant Integer InitialValues = 2 "Initial values for p, T or h, X";
    constant Integer SteadyState = 3 
      "Steady state (guess values for p, T or h, X)";
    constant Integer SteadyStateHydraulic = 4 
      "Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X";
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Modelica.Icons.TypeInteger(min=1, max=4);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.Init.NoInit 
            "NoInit (guess values for p, T or h, X)",
          choice=Modelica_Fluid.Types.Init.InitialValues 
            "InitialValues (initial values for p, T or h, X)",
          choice=Modelica_Fluid.Types.Init.SteadyState 
            "SteadyState (guess values for p, T or h, X)",
          choice=Modelica_Fluid.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (der(p)=0, guess value for p, initial values for T or h, X)"),
        Documentation(info="<html>
<p>
Integer type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.Init.</b></th><th><b>Meaning</b></th></tr>
<tr><td>NoInit (=1)</td>
    <td>No initial conditions (guess values for p, T or h, X)</td></tr>
 
<tr><td>InitialValues (=2)</td>
    <td>Initial values for p, T or h, X</td></tr>
 
<tr><td>SteadyState (=3)</td>
    <td>Steady state (guess values for p, T or h, X)</td></tr>
 
<tr><td>SteadyStateHydraulic (=4)</td>
    <td>Hydraulic steady state (der(p)=0), guess value for p, 
        initial values for T or h, X</td></tr>
</table>
</html>"));
      
    end Temp;
  end Init;
  
  package InitWithGlobalDefault 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    
    annotation (Documentation(info="<html>
 
</html>"));
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit = 1 
      "No initial conditions (guess values for p, T or h, X)";
    constant Integer InitialValues = 2 "Initial values for p, T or h, X";
    constant Integer SteadyState = 3 
      "Steady state (guess values for p, T or h, X)";
    constant Integer SteadyStateHydraulic = 4 
      "Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X";
    constant Integer UseGlobalFluidOption = 5 
      "Use initialization defined in fluidOptions component";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Modelica.Icons.TypeInteger(min=1, max=5);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.Init.NoInit 
            "NoInit (guess values for p, T or h, X)",
          choice=Modelica_Fluid.Types.Init.InitialValues 
            "InitialValues (initial values for p, T or h, X)",
          choice=Modelica_Fluid.Types.Init.SteadyState 
            "SteadyState (guess values for p, T or h, X)",
          choice=Modelica_Fluid.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (der(p)=0, guess value for p, initial values for T or h, X)",
          choice=Modelica_Fluid.Types.Init.UseGlobalFluidOption 
            "UseGlobalFluidOption (Use default defined in fluidOptions component)"),
        Documentation(info="<html>
<p>
Integer type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.InitWithGlobalDefault.</b></th><th><b>Meaning</b></th></tr>
<tr><td>NoInit (=1)</td>
    <td>No initial conditions (guess values for p, T or h, X)</td></tr>
 
<tr><td>InitialValues (=2)</td>
    <td>Initial values for p, T or h, X</td></tr>
 
<tr><td>SteadyState (=3)</td>
    <td>Steady state (guess values for p, T or h, X)</td></tr>
 
<tr><td>SteadyStateHydraulic (=4)</td>
    <td>Hydraulic steady state (der(p)=0), guess value for p, 
        initial values for T or h, X</td></tr>
 
<tr><td>UseGlobalFluidOption (=5)</td>
    <td>Use default initialization defined in \"inner\" 
        <a href=\"Modelica://Modelica_Fluid.Components.FluidOptions\">fluidOptions</a> 
        component</td></tr>
</table>
</html>"));
      
    end Temp;
  end InitWithGlobalDefault;
  
  package FlowDirection 
    "Type, constants and menu choices to define whether flow reversal is allowed, as temporary solution until enumerations are available" 
    
    annotation (Documentation(info="<html>
  
</html>"));
    extends Modelica.Icons.Enumeration;
    constant Integer Unidirectional = 1 
      "Fluid flows only from port_a to port_b";
    constant Integer Bidirectional = 2 
      "No restrictions on fluid flow (flow reversal possible)";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Modelica.Icons.TypeInteger(min=1, max=2);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.FlowDirection.Unidirectional 
            "Unidirectional (fluid flows only from port_a to port_b)",
          choice=Modelica_Fluid.Types.FlowDirection.Bidirectional 
            "Bidirectional (flow reversal possible)"),
        Documentation(info="<html>
<p>
Integer type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.FlowDirection.</b></th><th><b>Meaning</b></th></tr>
<tr><td>Unidirectional (=1)</td>
    <td>Fluid flows only from port_a to port_b.
        By this option, min and max values are set for
        port_a.m_flow and port_b.m_flow. This allows a
        Modelica translator to remove if-clauses of 
        the semiLinear(..) operator reducing the
        size of non-linear equation systems.</td></tr>
 
<tr><td>Bidirectional (=2)</td>
    <td>No restrictions on fluid flow (flow reversal possible)</td></tr>
 
</table>
</html>"));
    end Temp;
  end FlowDirection;
  
  package FlowDirectionWithGlobalDefault 
    "Type, constants and menu choices to define whether flow reversal is allowed, as temporary solution until enumerations are available" 
    
    extends Modelica.Icons.Enumeration;
    constant Integer Unidirectional = 1 "Fluid flows only in one direction";
    constant Integer Bidirectional = 2 
      "No restrictions on fluid flow (flow reversal possible)";
    constant Integer UseGlobalFluidOption = 3 
      "Use default flow direction defined in fluidOptions component";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Modelica.Icons.TypeInteger(min=1, max=3);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Unidirectional 
            "Unidirectional (fluid flows only from port_a to port_b)",
          choice=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional 
            "Bidirectional (flow reversal possible)",
          choice=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
            "UseGlobalFluidOption (use default defined in fluidOptions component)"),
        Documentation(info="<html>
<p>
Integer type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.FlowDirectionWithGlobalDefault.</b></th><th><b>Meaning</b></th></tr>
<tr><td>Unidirectional (=1)</td>
    <td>Fluid flows only from port_a to port_b.
        By this option, min and max values are set for
        port_a.m_flow and port_b.m_flow. This allows a
        Modelica translator to remove if-clauses of 
        the semiLinear(..) operator reducing the
        size of non-linear equation systems.</td></tr>
 
<tr><td>Bidirectional (=2)</td>
    <td>No restrictions on fluid flow (flow reversal possible)</td></tr>
 
<tr><td>UseGlobalFluidOption (=3)</td>
    <td>Use default flow direction defined in \"inner\"
<a href=\"Modelica://Modelica_Fluid.Components.FluidOptions\">fluidOptions</a> 
        component</td></tr>
 
</table>
</html>"));
      
    end Temp;
    annotation (Documentation(info="<html> 
  
</html>"));
  end FlowDirectionWithGlobalDefault;
  
  model CvTypes 
    "Type, constants and menu choices to define the choice of valve flow coefficient" 
    extends Modelica.Icons.Enumeration;
    annotation (preferedView="text");
    constant Integer Av = 0 "Av (metric) flow coefficient";
    constant Integer Kv = 1 "Kv (metric) flow coefficient";
    constant Integer Cv = 2 "Cv (US) flow coefficient";
    constant Integer OpPoint = 3 "Av defined by operating point";
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=0, max=3);
      annotation (Evaluate=true, choices(
        choice=Modelica_Fluid.Types.CvTypes.Av "Av (metric) flow coefficient",
        choice=Modelica_Fluid.Types.CvTypes.Kv "Kv (metric) flow coefficient",
        choice=Modelica_Fluid.Types.CvTypes.Cv "Cv (US) flow coefficient",
        choice=Modelica_Fluid.Types.CvTypes.OpPoint 
            "Av defined by nominal operating point"));
    end Temp;
    
  end CvTypes;
  
end Types;
