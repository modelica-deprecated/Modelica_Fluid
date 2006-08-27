package Sources 
  "Generic sources for fluid connectors to define fixed or prescribed ambient conditions" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  import SI = Modelica.SIunits;
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
            fillColor=69)), Text(extent=[-136, 144; 132, 82], string="%name")),
      Documentation(info="<html>
<p>
Model <b>FixedBoundary</b> defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure or boundary density.</li>
<li> Boundary temperature or boundary specific enthalpy.</li>
<li> Boundary mass fractions (only for multi-substance flow).</li>
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
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p 
      "Boundary pressure";
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T 
      "Boundary temperature";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X[Medium.nX](
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
            fillColor=69)), Text(extent=[-150,150; 150,100],   string="%name")),
      Documentation(info="<html>
<p>
Defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary temperature.</li>
<li> Boundary mass fractions (only for multi-substance flow).</li>
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
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p 
      "Boundary pressure";
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy h 
      "Boundary specific enthalpy";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X[
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
            fillColor=69)), Text(extent=[-136, 144; 132, 82], string="%name")),
      Documentation(info="<html>
<p>
Defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary specific enthalpy.</li>
<li> Boundary mass fractions (only for multi-substance flow).</li>
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
    parameter SI.Pressure p "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(p_in)==0)));
    parameter SI.Temperature T "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(T_in)==0)));
    parameter SI.MassFraction X[Medium.nX] = Medium.X_default 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(X_in)==0) or Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in(redeclare type SignalType = 
          SI.Pressure) "Prescribed boundary pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput T_in(
      redeclare type SignalType = SI.Temperature) 
      "Prescribed boundary temperature" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](redeclare type 
        SignalType = SI.MassFraction) "Prescribed boundary composition" 
      annotation (extent=[-140,-80; -100,-40]);
    annotation (defaultComponentName = "boundary_prescribed",
  Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=69,
            gradient=3,
            fillColor=69)), Text(extent=[-134,168; 134,106],  string="%name"),
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
<li> Prescribed boundary pressure via input signal <tt>p_in</tt>.</li>
<li> Prescribed boundary temperature via input signal <tt>T_in</tt>.</li>
<li> Prescribed boundary mass fractions via input signal <tt>X_in</tt> (only for multi-substance flow).</li>
</ul>
<p>If the connector are left unconnected, the corresponding prescribed values
are set by the parameters <tt>p</tt>, <tt>T</tt>, and <tt>X</tt>, respectively.
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
                                          Medium.singleState, true, X_in, "PrescribedBoundary_pTX");
    if cardinality(p_in)==0 then
      p_in = p;
    end if;
    if cardinality(T_in)==0 then
      T_in = T;
    end if;
    if cardinality(X_in)==0 then
      X_in = X;
    end if;
    medium.p = p_in;
    medium.T = T_in;
    medium.Xi = X_in[1:Medium.nXi];
  end PrescribedBoundary_pTX;
  
  model PrescribedBoundary_phX 
    "Boundary with prescribed pressure, specific enthalpy and composition" 
    extends Sources.BaseClasses.PartialSource;
    parameter SI.Pressure p "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(p_in)==0)));
    parameter SI.SpecificEnthalpy h "Fixed value of specific enthalpy" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(h_in)==0)));
    parameter SI.MassFraction X[Medium.nX] = Medium.X_default 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(X_in)==0) or Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in(
      redeclare type SignalType = SI.Pressure) "Prescribed boundary pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput h_in(
      redeclare type SignalType = SI.SpecificEnthalpy) 
      "Prescribed boundary specific enthalpy" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = SI.MassFraction) 
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
            fillColor=69)), Text(extent=[-134,168; 134,106],  string="%name"),
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
Defines values for boundary conditions:
</p>
<ul>
<li> Prescribed boundary pressure via input signal <tt>p_in</tt>.</li>
<li> Prescribed boundary specific enthalpy via input signal <tt>h_in</tt>.</li>
<li> Prescribed boundary mass fractions via input signal <tt>X_in</tt> (only for multi-substance flow).</li>
</ul>
<p>If the connector are left unconnected, the corresponding prescribed values
are set by the parameters <tt>p</tt>, <tt>h</tt>, and <tt>X</tt>, respectively.
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
                                          Medium.singleState, true, X_in, "PrescribedBoundary_phX");
    if cardinality(p_in)==0 then
      p_in = p;
    end if;
    if cardinality(h_in)==0 then
      h_in = h;
    end if;
    if cardinality(X_in)==0 then
      X_in = X;
    end if;
    medium.p = p_in;
    medium.h = h_in;
    medium.Xi = X_in[1:Medium.nXi];
  end PrescribedBoundary_phX;
  
  model PrescribedMassFlowRate_TX 
    "Ideal pump that produces a prescribed mass flow with prescribed temperature and mass fraction" 
    extends Sources.BaseClasses.PartialSource;
    parameter Medium.MassFlowRate m_flow = 0 
      "Fixed mass flow rate going out of the fluid port";
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T 
      "Fixed value of the fluid temperature";
    parameter Medium.MassFraction X[Medium.nX](quantity=Medium.substanceNames) = Medium.X_default 
      "Fixed value of the fluid composition" 
      annotation (Dialog(enable = Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in(
      redeclare type SignalType = SI.MassFlowRate) "Prescribed mass flow rate" 
      annotation (extent=[-113,40; -73,80]);
    Modelica.Blocks.Interfaces.RealInput T_in(
      redeclare type SignalType = SI.Temperature) 
      "Prescribed fluid temperature" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = SI.MassFraction) 
      "Prescribed fluid composition" 
      annotation (extent=[-112,-81; -72,-41]);
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
        Text(extent=[-150,160; 150,110],  string="%name"),
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
<li> Prescribed mass flow rate via input signal <tt>m_flow_in</tt>.</li>
<li> Prescribed temperature via input signal <tt>T_in</tt>.</li>
<li> Prescribed mass fractions via input signal <tt>X_in</tt> (only for multi-substance flow) .</li>
</ul>
<p>If the connector are left unconnected, the corresponding prescribed values
are set by the parameters <tt>m_flow</tt>, <tt>T</tt>, and <tt>X</tt>, respectively.
<p>
Note, that temperature
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
  outer Modelica_Fluid.Ambient ambient "Ambient conditions";
  equation 
    Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                           Medium.singleState, true, X, "PrescribedMassFlowRate_TX");
    if cardinality(m_flow_in)==0 then
      m_flow_in = m_flow;
    end if;
    if cardinality(T_in)==0 then
      T_in = T;
    end if;
    if cardinality(X_in)==0 then
      X_in = X;
    end if;
    port.m_flow = -m_flow_in;
    medium.T = T_in;
    medium.Xi = X_in[1:Medium.nXi];
  end PrescribedMassFlowRate_TX;
  
  model PrescribedMassFlowRate_hX 
    "Ideal pump that produces a prescribed mass flow with prescribed specific enthalpy and mass fraction" 
    extends Sources.BaseClasses.PartialSource;
    parameter Medium.MassFlowRate m_flow = 0 
      "Fixed mass flow rate going out of the fluid port";
    parameter Medium.SpecificEnthalpy h 
      "Fixed value of the fluid specific enthalpy";
    parameter Medium.MassFraction X[Medium.nX](quantity=Medium.substanceNames) = Medium.X_default 
      "Fixed value of the fluid composition" 
      annotation (Dialog(enable=Medium.nXi>0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in(
      redeclare type SignalType = SI.MassFlowRate) "Prescribed mass flow rate" 
      annotation (extent=[-113,40; -73,80]);
    Modelica.Blocks.Interfaces.RealInput h_in(
      redeclare type SignalType = SI.SpecificEnthalpy) 
      "Prescribed fluid specific enthalpy" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = SI.MassFraction) 
      "Prescribed fluid composition" 
      annotation (extent=[-113,-80; -73,-40]);
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
        Text(extent=[-150,160; 150,110],  string="%name")),
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
<li> Prescribed mass flow rate via input signal <tt>m_flow_in</tt>.</li>
<li> Prescribed specific enthalpy via input signal <tt>h_in</tt>.</li>
<li> Prescribed mass fractions via input signal <tt>X_in</tt> (only for multi-substance flow) .</li>
</ul>
<p>If the connector are left unconnected, the corresponding prescribed values
are set by the parameters <tt>m_flow</tt>, <tt>h</tt>, and <tt>X</tt>, respectively.
<p>
Note, that specific enthalpy
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
  equation 
    Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                           Medium.singleState, true, X, "PrescribedMassFlowRate_hX");
    if cardinality(m_flow_in)==0 then
      m_flow_in = m_flow;
    end if;
    if cardinality(h_in)==0 then
      h_in = h;
    end if;
    if cardinality(X_in)==0 then
      X_in = X;
    end if;
    port.m_flow = -m_flow_in;
    medium.h = h_in;
    medium.Xi = X_in[1:Medium.nXi];
  end PrescribedMassFlowRate_hX;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
  partial model PartialSource 
      "Partial component source with one fluid connector (default" 
      import Modelica.Constants;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
       annotation (choicesAllMatching=true);
    parameter Types.SourceFlowDirection.Temp flowDirection=
                     Types.SourceFlowDirection.Bidirectional 
        "Uni- or bidirectional flow component" 
                                             annotation(Evaluate=true, Dialog(tab="Advanced"));
      
    Medium.BaseProperties medium "Medium in the source";
    Interfaces.FluidPort_b port(redeclare package Medium = Medium,
                     m_flow(max=if flowDirection==Types.SourceFlowDirection.OutOfPort then 0 else 
                                     +Constants.inf,
                            min=if flowDirection==Types.SourceFlowDirection.InToPort then 0 else 
                                     -Constants.inf)) 
      annotation (extent=[90,-10; 110,10],    rotation=0);
      
  equation 
    port.p = medium.p;
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
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
  end PartialSource;
  end BaseClasses;
end Sources;
