package Types 
  package Init 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    
    annotation (preferedView="text");
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
      extends Integer(min=1, max=4);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.Init.NoInit 
            "NoInit (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.InitialValues 
            "InitialValues (initial values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyState 
            "SteadyState (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (der(p)=0, guess value for p, initial values for T or h, X)"));
    end Temp;
  end Init;
  
  package InitWithGlobalDefault 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit = 1 
      "No initial conditions (guess values for p, T or h, X)";
    constant Integer InitialValues = 2 "Initial values for p, T or h, X";
    constant Integer SteadyState = 3 
      "Steady state (guess values for p, T or h, X)";
    constant Integer SteadyStateHydraulic = 4 
      "Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X";
    constant Integer UseEnvironmentOption = 5 
      "Use initialization defined in environment component";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=5);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.Init.NoInit 
            "NoInit (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.InitialValues 
            "InitialValues (initial values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyState 
            "SteadyState (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (der(p)=0, guess value for p, initial values for T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.UseEnvironmentOption 
            "UseEnvironmentOption (use initialization defined in environment component)"));
    end Temp;
  end InitWithGlobalDefault;
  
  package Flow 
    "Type, constants and menu choices to define whether flow reversal is allowed, as temporary solution until enumerations are available" 
    
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer Unidirectional = 1 "Fluid flows only in one direction";
    constant Integer Bidirectional = 2 
      "No restrictions on fluid flow (flow reversal possible)";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=2);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Unidirectional 
            "Unidirectional (fluid flows only in one direction)",
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Bidirectional 
            "Bidirectional (flow reversal possible)"));
    end Temp;
  end Flow;
  
  package FlowWithGlobalDefault 
    "Type, constants and menu choices to define whether flow reversal is allowed, as temporary solution until enumerations are available" 
    
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer Unidirectional = 1 "Fluid flows only in one direction";
    constant Integer Bidirectional = 2 
      "No restrictions on fluid flow (flow reversal possible)";
    constant Integer UseEnvironmentOption = 3 
      "Use FlowReversal defined in environment component";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=3);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Unidirectional 
            "Unidirectional (fluid flows only in one direction)",
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Bidirectional 
            "Bidirectional (flow reversal possible)",
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.UseEnvironmentOption 
            "UseEnvironmentOption (use FlowReversal defined in environment component)"));
    end Temp;
  end FlowWithGlobalDefault;
  
  package OrificeCharacteristics "Functions for valve characteristics" 
    
    annotation (preferedView="info",
  Documentation(info="<html>
<p>
This package provides functions that compute a 
<b>constant pressure loss factor</b> for 
<b>turbulent flow</b> based solely on geometrical
information of the corresponding component.
All functions have to be derived from function 
<b>partialTurbulentLossFactor</b>.
</p>
<p>
The turbulent pressure loss factor is defined
according to the following equations:
</p>
<pre>
  dp       = port_a.p - port_b.p     \"pressure drop from port_a to port_b\"
  m_flow_a = port_a.m_flow           \"mass flow rate into port_a\"
  d_a      = d_a(port_a.p, port_a.h) \"density at port_a\" 
  v_a      = m_flow_a/(d_a*A_a)      \"velocity at port_a\"
&nbsp;
  // It is assumed that m_flow_a > 0
  dp = 0.5*d_a*zeta_a*v_a^2
     = 0.5*d_a*zeta_a*m_flow_a^2/(d_a*A_a)^2
     = (0.5*zeta_a/A_a^2) * (1/d_a) * m_flow_a^2
</pre>
<p>
When the flow is reversed, the pressure loss equation is defined
in the other way round, using the quantities at port_b.
This gives the following overall description:
</p>
<pre>
  dp = <b>if</b> m_flow_a &gt; 0 <b>then</b>
           k_a * 1/d_a * m_flow_a^2
       <b>else</b>
          -k_b * 1/d_b * m_flow_b^2
&nbsp;
       k_a = 0.5*zeta_a/A_a^2
       k_b = 0.5*zeta_b/A_b^2
</pre>
<p>
This functions returns the pressure loss factors k_a and k_b
as a vector with two elements:
</p>
<pre>
  k[1]: loss factor k_a if m_flow_a &gt; 0
  k[2]: loss factor k_b if m_flow_a &lt; 0
</pre>
<p>
If a pressure loss factor is only known for one flow
direction, m_flow_a &gt; 0, the component can only be
used, if the flow is most of the time in this direction.
One might define k[2]=k[1], if only for short time periods
also a reversal flow appears and then neglect the error
for the flow reversal and at small mass flow rates
when the flow is laminar.
</p>
</html>"));
    
    model turbulentLossFactor "Constant loss factors for turbulent flow" 
      parameter Real k_a "Loss factor if fluid flows from port_a to port_b";
      parameter Real k_b=k_a "Loss factor if fluid flows from port_b to port_a";
    end turbulentLossFactor;
    
    model suddenExpansion "Suddenly expanding area" 
      import SI = Modelica.SIunits;
      parameter SI.Area A_a "Area at port_a";
      parameter SI.Area A_b "Area at port_b (A_b > A_a required)";
      extends turbulentLossFactor(
              final k_a = (1 - A_a/A_b)^2/A_a^2,
              final k_b = 0.5*(1 - A_a/A_b)^0.75/A_a^2);
      annotation (Icon(
          Rectangle(extent=[-100,40; 0,-40], style(
              color=7,
              rgbcolor={255,255,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Rectangle(extent=[0,100; 100,-100], style(
              color=7,
              rgbcolor={255,255,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                100; 0,40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-76,4; 34,4], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Polygon(points=[34,16; 34,-10; 74,4; 34,16], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1))), Diagram(
          Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                100; 0,40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Rectangle(extent=[-100,40; 0,-40], style(
              color=7,
              rgbcolor={255,255,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Rectangle(extent=[0,100; 100,-100], style(
              color=7,
              rgbcolor={255,255,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                100; 0,40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-60,-40; -60,40], style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1)),
          Text(
            extent=[-50,16; -26,-10],
            style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1),
            string="A_a"),
          Line(points=[34,-100; 34,100], style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1)),
          Text(
            extent=[54,16; 78,-10],
            style(
              color=3,
              rgbcolor={0,0,255},
              arrow=3,
              fillColor=3,
              rgbfillColor={0,0,255},
              fillPattern=1),
            string="A_b")));
    equation 
      
    end suddenExpansion;
    
  end OrificeCharacteristics;
end Types;
