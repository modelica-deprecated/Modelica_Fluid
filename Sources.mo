within Modelica_Fluid;
package Sources 
  "Generic sources for fluid connectors to define fixed or prescribed ambient conditions" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  model FixedBoundary "Boundary source component" 
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean use_p=true "select p or d" 
      annotation (Evaluate = true,
                  Dialog(group = "Boundary pressure or Boundary density"));
    parameter Medium.AbsolutePressure p "Boundary pressure" 
      annotation (Dialog(group = "Boundary pressure or Boundary density",
                         enable = use_p));
    parameter Medium.Density d=1000 "Boundary density" 
      annotation (Dialog(group = "Boundary pressure or Boundary density",
                         enable=not use_p));
    parameter Boolean use_T=true "select T or h" 
      annotation (Evaluate = true,
                  Dialog(group = "Boundary temperature or Boundary specific enthalpy"));
    parameter Medium.Temperature T "Boundary temperature" 
      annotation (Dialog(group = "Boundary temperature or Boundary specific enthalpy",
                         enable = use_T));
    parameter Medium.SpecificEnthalpy h "Boundary specific enthalpy" 
      annotation (Dialog(group="Boundary temperature or Boundary specific enthalpy",
                  enable = not use_T));
    parameter Medium.MassFraction X[Medium.nX](
         quantity=Medium.substanceNames)=Medium.X_default 
      "Boundary mass fractions m_i/m" 
      annotation (Dialog(group = "Only for multi-substance flow", enable=Medium.nXi > 0));
    
    annotation (defaultComponentName = "Boundary_fixed",
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-150,110; 150,150], string="%name")),
      Documentation(info="<html>
<p>
Model <b>FixedBoundary</b> defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure or boundary density.</li>
<li> Boundary temperature or boundary specific enthalpy.</li>
<li> Boundary composition (only for multi-substance flow).</li>
</ul>
<p>
Note, that boundary temperature, density, specific enthalpy
and mass fractions have only an effect if the mass flow
is from the Boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
    
  equation 
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                          Medium.singleState, use_p, X,
                                          "FixedBoundary");
    if use_p or Medium.singleState then
      medium.p = p;
    else
      medium.d = d;
    end if;
    if use_T then
      medium.T = T;
    else
      medium.h = h;
    end if;
    
    medium.Xi = X[1:Medium.nXi];
  end FixedBoundary;
  
  model FixedBoundary_pTX 
    "Boundary pressure, temperature and mass fraction source" 
    extends Sources.BaseClasses.PartialSource;
    parameter Medium.AbsolutePressure p "Boundary pressure";
    parameter Medium.Temperature T "Boundary temperature";
    parameter Medium.MassFraction X[Medium.nX](
         quantity=Medium.substanceNames) = Medium.X_default 
      "Boundary mass fractions m_i/m" 
      annotation (Dialog(group = "Only for multi-substance flow",
                  enable=Medium.nXi > 0));
    annotation (defaultComponentName = "boundary_fixed",
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-150,110; 150,150],   string="%name")),
      Documentation(info="<html>
<p>
Defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary temperature.</li>
<li> Boundary composition (only for multi-substance flow).</li>
</ul>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
  equation 
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                          Medium.singleState, true, X, "FixedBoundary_pTX");
    medium.p   = p;
    medium.T   = T;
    medium.Xi = X[1:Medium.nXi];
  end FixedBoundary_pTX;
  
  model FixedBoundary_phX 
    "Boundary pressure, specific enthalpy and mass fraction source" 
    extends Sources.BaseClasses.PartialSource;
    parameter Medium.AbsolutePressure p "Boundary pressure";
    parameter Medium.SpecificEnthalpy h "Boundary specific enthalpy";
    parameter Medium.MassFraction X[
      Medium.nX](quantity=Medium.substanceNames) = Medium.X_default 
      "Boundary mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nXi > 0));
    annotation (defaultComponentName = "boundary_fixed",
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-150,110; 150,150], string="%name")),
      Documentation(info="<html>
<p>
Defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary specific enthalpy.</li>
<li> Boundary composition (only for multi-substance flow).</li>
</ul>
<p>
Note, that boundary specific enthalpy
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
  equation 
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                          Medium.singleState, true, X, "FixedBoundary_phX");
    medium.p = p;
    medium.h = h;
    medium.Xi = X[1:Medium.nXi];
  end FixedBoundary_phX;
  
  annotation (Documentation(info="<html>
<p>
Package <b>Sources</b> contains generic sources for fluid connectors
to define fixed or prescribed ambient conditions.
</p>
</html>"));
  model PrescribedBoundary_pTX 
    "Boundary with prescribed pressure, temperature and composition" 
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean usePressureInput = false 
      "Get the pressure from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useTemperatureInput= false 
      "Get the temperature from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false 
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.AbsolutePressure p = Medium.reference_p 
      "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = not usePressureInput));
    parameter Medium.Temperature T = Medium.reference_T 
      "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = not useTemperatureInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in(
      redeclare type SignalType = Medium.AbsolutePressure) if usePressureInput 
      "Prescribed boundary pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput T_in(
      redeclare type SignalType = Medium.Temperature) if useTemperatureInput 
      "Prescribed boundary temperature" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) if useCompositionInput 
      "Prescribed boundary composition" 
      annotation (extent=[-140,-80; -100,-40]);
  protected 
    Modelica.Blocks.Interfaces.RealInput p_in_internal(
      redeclare type SignalType = Medium.AbsolutePressure) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput T_in_internal(
      redeclare type SignalType = Medium.Temperature) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](
     redeclare type SignalType = Medium.MassFraction) 
      "Needed to connect to conditional connector";
    annotation (defaultComponentName = "boundary_prescribed",
  Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-150,110; 150,150],  string="%name"),
        Line(points=[-100,60; -80,60], style(color=3, rgbcolor={0,0,255})),
        Line(points=[-100,-60; -80,-60], style(color=3, rgbcolor={0,0,255})),
        Text(
          extent=[-146,110; -62,70],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="p"),
        Text(
          extent=[-160,-22; -58,-62],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="X"),
        Text(
          extent=[-158,44; -56,4],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="T")),
      Documentation(info="<html>
<p>
Defines prescribed values for boundary conditions:
</p>
<ul>
<li> Prescribed boundary pressure.</li>
<li> Prescribed boundary temperature.</li>
<li> Prescribed boundary composition (only for multi-substance flow).</li>
</ul>
<p>If <tt>usePressureInput</tt> is false (default option), the <tt>p</tt> parameter
is used as boundary pressure, and the <tt>p_in</tt> input connector is disabled; if <tt>usePressureInput</tt> is true, then the <tt>p</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"),
      Diagram);
  equation 
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedBoundary_pTX");
    connect(p_in, p_in_internal);
    connect(T_in, T_in_internal);
    connect(X_in, X_in_internal);
    if not usePressureInput then
      p_in_internal = p;
    end if;
    if not useTemperatureInput then
      T_in_internal = T;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    medium.p = p_in_internal;
    medium.T = T_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedBoundary_pTX;
  
  model PrescribedBoundary_phX 
    "Boundary with prescribed pressure, specific enthalpy and composition" 
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean usePressureInput = false 
      "Get the pressure from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useEnthalpyInput= false 
      "Get the specific enthalpy from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false 
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.AbsolutePressure p = Medium.reference_p 
      "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = not usePressureInput));
    parameter Medium.SpecificEnthalpy h = Medium.h_default 
      "Fixed value of specific enthalpy" 
      annotation (Evaluate = true,
                  Dialog(enable = not useEnthalpyInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in(
      redeclare type SignalType = Medium.AbsolutePressure) if usePressureInput 
      "Prescribed boundary pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput h_in(
      redeclare type SignalType = Medium.SpecificEnthalpy) if useEnthalpyInput 
      "Prescribed boundary specific enthalpy" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) if useCompositionInput 
      "Prescribed boundary composition" 
      annotation (extent=[-140,-80; -100,-40]);
    annotation (defaultComponentName = "boundary_prescribed",
  Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-150,110; 150,150],  string="%name"),
        Line(points=[-100,60; -66,60], style(color=3, rgbcolor={0,0,255})),
        Line(points=[-100,-60; -66,-60], style(color=3, rgbcolor={0,0,255})),
        Text(
          extent=[-148,120; -70,80],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="p"),
        Text(
          extent=[-152,-86; -50,-126],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="h")),
      Documentation(info="<html>
<p>
Defines prescribed values for boundary conditions:
</p>
<ul>
<li> Prescribed boundary pressure.</li>
<li> Prescribed boundary temperature.</li>
<li> Prescribed boundary composition (only for multi-substance flow).</li>
</ul>
<p>If <tt>usePressureInput</tt> is false (default option), the <tt>p</tt> parameter
is used as boundary pressure, and the <tt>p_in</tt> input connector is disabled; if <tt>usePressureInput</tt> is true, then the <tt>p</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the specific enthalpy and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
  protected 
    Modelica.Blocks.Interfaces.RealInput p_in_internal(
      redeclare type SignalType = Medium.AbsolutePressure) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput h_in_internal(
      redeclare type SignalType = Medium.SpecificEnthalpy) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) 
      "Needed to connect to conditional connector";
  equation 
    Modelica_Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedBoundary_phX");
    connect(p_in, p_in_internal);
    connect(h_in, h_in_internal);
    connect(X_in, X_in_internal);
    if not usePressureInput then
      p_in_internal = p;
    end if;
    if not useEnthalpyInput then
      h_in_internal = h;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    medium.p = p_in_internal;
    medium.h = h_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedBoundary_phX;
  
  model PrescribedMassFlowRate_TX 
    "Ideal flow source that produces a prescribed mass flow with prescribed temperature and mass fraction" 
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean useFlowRateInput = false 
      "Get the mass flow rate from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useTemperatureInput= false 
      "Get the temperature from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false 
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.MassFlowRate m_flow = 0 
      "Fixed mass flow rate going out of the fluid port" 
      annotation (Evaluate = true,
                  Dialog(enable = not useFlowRateInput));
    parameter Medium.Temperature T = Medium.reference_T 
      "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = not useTemperatureInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in(
      redeclare type SignalType = Medium.MassFlowRate) if useFlowRateInput 
      "Prescribed mass flow rate" 
      annotation (extent=[-113,40; -73,80]);
    Modelica.Blocks.Interfaces.RealInput T_in(
      redeclare type SignalType = Medium.Temperature) if useTemperatureInput 
      "Prescribed fluid temperature" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) if useCompositionInput 
      "Prescribed fluid composition" 
      annotation (extent=[-112,-81; -72,-41]);
  protected 
    Modelica.Blocks.Interfaces.RealInput m_flow_in_internal(
      redeclare type SignalType = Medium.MassFlowRate) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput T_in_internal(
      redeclare type SignalType = Medium.Temperature) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) 
      "Needed to connect to conditional connector";
    annotation (defaultComponentName = "massFlowRate",
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1,1],
        component=[20, 20],
        scale=0),
      Icon(
        Rectangle(extent=[20, 60; 100, -60], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[38, 40; 100, -40], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
        Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
              fillColor=73)),
        Text(
          extent=[-54, 32; 16, -30],
          style(color=41, fillColor=41),
          string="m"),
        Text(extent=[-150,110; 150,150],  string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Text(
          extent=[-194,112; -54,80],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="m_flow"),
        Text(
          extent=[-100,14; -60,-20],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="T"),
        Text(
          extent=[-144,-90; -24,-118],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="X")),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram,
      Documentation(info="<html>
<p>
Models an ideal flow source, with prescribed values of flow rate, temperature and composition:
</p>
<ul>
<li> Prescribed mass flow rate.</li>
<li> Prescribed temperature.</li>
<li> Prescribed composition (only for multi-substance flow) .</li>
</ul>
<p>If <tt>useFlowRateInput</tt> is false (default option), the <tt>m_flow</tt> parameter
is used as boundary pressure, and the <tt>m_flow_in</tt> input connector is disabled; if <tt>useFlowRateInput</tt> is true, then the <tt>m_flow</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary flow rate, do not have an effect.
</p>
</html>"));
  equation 
    Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedMassFlowRate_TX");
    connect(m_flow_in, m_flow_in_internal);
    connect(T_in, T_in_internal);
    connect(X_in, X_in_internal);
    if not useFlowRateInput then
      m_flow_in_internal = m_flow;
    end if;
    if not useTemperatureInput then
      T_in_internal = T;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    port.m_flow = -m_flow_in_internal;
    medium.T = T_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedMassFlowRate_TX;
  
  model PrescribedMassFlowRate_hX 
    "Ideal flow source that produces a prescribed mass flow with prescribed specific enthalpy and mass fraction" 
    extends Sources.BaseClasses.PartialSource;
    parameter Boolean useFlowRateInput = false 
      "Get the mass flow rate from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useEnthalpyInput= false 
      "Get the specific enthalpy from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Boolean useCompositionInput = false 
      "Get the composition from the input connector" 
      annotation(Evaluate=true, HideResult=true, choices(__Dymola_checkBox=true));
    parameter Medium.MassFlowRate m_flow = 0 
      "Fixed mass flow rate going out of the fluid port" 
      annotation (Evaluate = true,
                  Dialog(enable = not useFlowRateInput));
    parameter Medium.SpecificEnthalpy h = Medium.h_default 
      "Fixed value of specific enthalpy" 
      annotation (Evaluate = true,
                  Dialog(enable = not useEnthalpyInput));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (not useCompositionInput) and Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in(
      redeclare type SignalType = Medium.MassFlowRate) if useFlowRateInput 
      "Prescribed mass flow rate" 
      annotation (extent=[-113,40; -73,80]);
    Modelica.Blocks.Interfaces.RealInput h_in(
      redeclare type SignalType = Medium.SpecificEnthalpy) if useEnthalpyInput 
      "Prescribed fluid specific enthalpy" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) if useCompositionInput 
      "Prescribed fluid composition" 
      annotation (extent=[-113,-80; -73,-40]);
  protected 
    Modelica.Blocks.Interfaces.RealInput m_flow_in_internal(
      redeclare type SignalType = Medium.MassFlowRate) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput h_in_internal(
      redeclare type SignalType = Medium.SpecificEnthalpy) 
      "Needed to connect to conditional connector";
    Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](
      redeclare type SignalType = Medium.MassFraction) 
      "Needed to connect to conditional connector";
    annotation (defaultComponentName = "massFlowRate",
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1,1],
        component=[20, 20],
        scale=0),
      Icon(
        Rectangle(extent=[20, 60; 100, -60], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[38, 40; 100, -40], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
        Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
              fillColor=73)),
        Text(
          extent=[-54, 32; 16, -30],
          style(color=41, fillColor=41),
          string="m"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Text(
          extent=[-194,115; -54,83],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="m_flow"),
        Text(
          extent=[-100,15; -60,-19],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="T"),
        Text(
          extent=[-145,-85; -25,-113],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="X"),
        Text(extent=[-150,110; 150,150],  string="%name")),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram,
      Documentation(info="<html>
<p>
Models an ideal flow source, with prescribed values of flow rate, temperature and composition:
</p>
<ul>
<li> Prescribed mass flow rate.</li>
<li> Prescribed specific enthalpy.</li>
<li> Prescribed composition (only for multi-substance flow) .</li>
</ul>
<p>If <tt>useFlowRateInput</tt> is false (default option), the <tt>m_flow</tt> parameter
is used as boundary pressure, and the <tt>m_flow_in</tt> input connector is disabled; if <tt>useFlowRateInput</tt> is true, then the <tt>m_flow</tt> parameter is ignored, and the value provided by the input connector is used instead.</p> 
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature
and mass fractions have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary flow rate, do not have an effect.
</p>
</html>"));
  equation 
    Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
      Medium.singleState, true, X_in_internal, "PrescribedMassFlowRate_hX");
    connect(m_flow_in, m_flow_in_internal);
    connect(h_in, h_in_internal);
    connect(X_in, X_in_internal);
    if not useFlowRateInput then
      m_flow_in_internal = m_flow;
    end if;
    if not useEnthalpyInput then
      h_in_internal = h;
    end if;
    if not useCompositionInput then
      X_in_internal = X;
    end if;
    port.m_flow = -m_flow_in_internal;
    medium.h = h_in_internal;
    medium.Xi = X_in_internal[1:Medium.nXi];
  end PrescribedMassFlowRate_hX;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
  partial model PartialSource 
      "Partial component source with one fluid connector" 
      import Modelica.Constants;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
       annotation (choicesAllMatching=true);
    parameter Types.SourceFlowDirection.Temp flowDirection=
                     Types.SourceFlowDirection.Bidirectional 
        "Allowed flow direction"             annotation(Evaluate=true, Dialog(tab="Advanced"));
      
    Medium.BaseProperties medium "Medium in the source";
    Modelica_Fluid.Interfaces.FluidPort_b port(
                                redeclare package Medium = Medium,
                     m_flow(max=if flowDirection==Types.SourceFlowDirection.OutOfPort then 0 else 
                                     +Constants.inf,
                            min=if flowDirection==Types.SourceFlowDirection.InToPort then 0 else 
                                     -Constants.inf)) 
      annotation (extent=[90,-10; 110,10],    rotation=0);
    annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_flow) and the mass flow rates of the
     substances (= port.mX_flow) depend on the direction of the mass flow rate.</li>
</ul>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  equation 
    port.p = medium.p;
    port.h_outflow  = medium.h;
    port.Xi_outflow = medium.Xi;
  end PartialSource;
  end BaseClasses;
end Sources;
