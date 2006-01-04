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
  
  package initOptions 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
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
  end initOptions;
  
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
  
  package ValveCharacteristics "Functions for valve characteristics" 
    partial function baseFun "Base class for valve characteristics" 
      extends Modelica.Icons.Function;
      input Real pos "Stem position (per unit)";
      output Real rc "Relative coefficient (per unit)";
    end baseFun;
    
    function linear "Linear characteristic" 
      extends baseFun;
    algorithm 
      rc := pos;
    end linear;
    
    function one "Constant characteristic" 
      extends baseFun;
    algorithm 
      rc := 1;
    end one;
    
    function quadratic "Quadratic characteristic" 
      extends baseFun;
    algorithm 
      rc := pos*pos;
    end quadratic;
    
    function equalPercentage "Equal percentage characteristic" 
      extends baseFun;
      input Real rangeability = 20 "Rangeability";
      input Real delta = 0.01;
    algorithm 
      rc := if pos > delta then rangeability^(pos-1) else 
              pos/delta*rangeability^(delta-1);
      annotation (Documentation(info="<html>
This characteristic is such that the relative change of the flow coefficient is proportional to the change in the stem position:
<p> d(rc)/d(pos) = k d(pos).
<p> The constant k is expressed in terms of the rangeability, i.e. the ratio between the maximum and the minimum useful flow coefficient:
<p> rangeability = exp(k) = rc(1.0)/rc(0.0).
<p> The theoretical characteristic has a non-zero opening when pos = 0; the implemented characteristic is modified so that the valve closes linearly when pos &lt delta.
</html>"));
    end equalPercentage;
    
  end ValveCharacteristics;
  
  package PumpCharacteristics "Functions for pump characteristics" 
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
    
    partial function baseFlow "Base class for pump flow characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output SI.Height head "Pump head";
    end baseFlow;
    
    partial function basePower 
      "Base class for pump power consumption characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output SI.Power consumption "Power consumption";
    end basePower;
    
    partial function baseEfficiency "Base class for efficiency characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output Real eta "Efficiency";
    end baseEfficiency;
    
    function linearFlow "Linear flow characteristic" 
      extends baseFlow;
      input SI.VolumeFlowRate q_nom[2] 
        "Volume flow rate for two operating points (single pump)";
      input SI.Height head_nom[2] "Pump head for two operating points";
    protected 
      constant Real g = Modelica.Constants.g_n;
      /* Linear system to determine the coefficients:
  head_nom[1]*g = c[1] + q_nom[1]*c[2];
  head_nom[2]*g = c[1] + q_nom[2]*c[2];
  */
      Real c[2] = Modelica.Math.Matrices.solve([ones(2),q_nom],head_nom*g) 
        "Coefficients of linear head curve";
    algorithm 
      // Flow equation: head * g = q*c[1] + c[2];
      head := 1/g * (c[1] + q_flow*c[2]);
    end linearFlow;
    
    function quadraticFlow "Quadratic flow characteristic" 
      extends baseFlow;
      input SI.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)";
      input SI.Height head_nom[3] "Pump head for three operating points";
    protected 
      constant Real g = Modelica.Constants.g_n;
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2} 
        "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  head_nom[1]*g = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  head_nom[2]*g = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  head_nom[3]*g = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3), q_nom, q_nom2],head_nom*g) 
        "Coefficients of quadratic head curve";
    algorithm 
      // Flow equation: head * g = c[1] + q_flow*c[2] + q_flow^2*c[3];
      head := 1/g * (c[1] + q_flow*c[2] + q_flow^2*c[3]);
    end quadraticFlow;
    
    function polynomialFlow "Polynomial flow characteristic" 
      extends baseFlow;
      input SI.VolumeFlowRate q_nom[:] 
        "Volume flow rate for N operating points (single pump)";
      input SI.Height head_nom[:] "Pump head for N operating points";
    protected 
      constant Real g = Modelica.Constants.g_n;
      Integer N = size(q_nom,1) "Number of nominal operating points";
      Real q_nom_pow[N,N] = {{q_nom[j]^(i-1) for j in 1:N} for i in 1:N} 
        "Rows: different operating points; columns: increasing powers";
      /* Linear system to determine the coefficients (example N=3):
  head_nom[1]*g = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  head_nom[2]*g = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  head_nom[3]*g = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[N] = Modelica.Math.Matrices.solve(q_nom_pow,head_nom*g) 
        "Coefficients of polynomial head curve";
    algorithm 
      // Flow equation (example N=3): head * g = c[1] + q_flow*c[2] + q_flow^2*c[3];
      // Note: the implementation is numerically efficient only for low values of Na
      head := 1/g * sum(q_flow^(i-1)*c[i] for i in 1:N);
    end polynomialFlow;
    
    function constantEfficiency "Constant efficiency characteristic" 
       extends baseEfficiency;
       input Real eta_nom "Nominal efficiency";
    algorithm 
      eta := eta_nom;
    end constantEfficiency;
    
    function quadraticPower "Quadratic power consumption characteristic" 
      extends basePower;
      input SI.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)";
      input SI.Power W_nom[3] "Power consumption for three operating points";
    protected 
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2} 
        "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  W_nom[1]*g = c[1] + q_nom[1]*c[2] + q_nom[1]^2*c[3];
  W_nom[2]*g = c[1] + q_nom[2]*c[2] + q_nom[2]^2*c[3];
  W_nom[3]*g = c[1] + q_nom[3]*c[2] + q_nom[3]^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([ones(3),q_nom,q_nom2],W_nom) 
        "Coefficients of quadratic power consumption curve";
    algorithm 
      consumption := c[1] + q_flow*c[2] + q_flow^2*c[3];
    end quadraticPower;
    
  end PumpCharacteristics;
end Types;
