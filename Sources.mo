package Sources "Generic fluid sources" 
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  model FixedAmbient "Ambient source component" 
    extends Interfaces.PartialSource;
    parameter Boolean use_p=true "select p or d" 
      annotation (Evaluate = true,
                  Dialog(group = "Ambient pressure or ambient density"));
    parameter Medium.AbsolutePressure p = Medium.reference_p "Ambient pressure"
      annotation (Dialog(group = "Ambient pressure or ambient density",
                         enable = use_p));
    parameter Medium.Density d=1000 "Ambient density" 
      annotation (Dialog(group = "Ambient pressure or ambient density",
                         enable=not use_p));
    parameter Boolean use_T=true "select T or h" 
      annotation (Evaluate = true,
                  Dialog(group = "Ambient temperature or ambient specific enthalpy"));
    parameter Medium.Temperature T = Modelica.SIunits.Conversions.from_degC(20) 
      "Ambient temperature" 
      annotation (Dialog(group = "Ambient temperature or ambient specific enthalpy",
                         enable = use_T));
    parameter Medium.SpecificEnthalpy h = 1.e4 "Ambient specific enthalpy" 
      annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                  enable = not use_T));
    parameter Medium.MassFraction X[Medium.nX](
         quantity=Medium.substanceNames)=Medium.reference_X 
      "Ambient mass fractions m_i/m" 
      annotation (Dialog(group = "Only for multi-substance flow", enable=Medium.nXi > 0));
    
    annotation (
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
Model <b>FixedAmbient</b> defines constant values for ambient conditions:
</p>
<ul>
<li> Ambient pressure or ambient density.</li>
<li> Ambient temperature or ambient specific enthalpy.</li>
<li> Ambient mass fractions (only for multi-substance flow).</li>
</ul>
<p>
Note, that ambient temperature, density, specific enthalpy
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
    
  equation 
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      use_p, X);
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
  end FixedAmbient;
  
  model FixedAmbient_pTX 
    "Ambient pressure, temperature and mass fraction source" 
    extends Interfaces.PartialSource;
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p=
        Medium.reference_p "Ambient pressure";
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T=
        Modelica.SIunits.Conversions.from_degC(20) "Ambient temperature";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X[Medium.nX](
         quantity=Medium.substanceNames) = Medium.reference_X 
      "Ambient mass fractions m_i/m" 
      annotation (Dialog(group = "Only for multi-substance flow",
                  enable=Medium.nXi > 0));
    annotation (
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
Model <b>FixedAmbient_pt</b> defines constant values for ambient conditions:
</p>
<ul>
<li> Ambient pressure.</li>
<li> Ambient temperature.</li>
<li> Ambient mass fractions (only for multi-substance flow).</li>
</ul>
<p>
Note, that ambient temperature
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
  equation 
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      true, X);
    medium.p   = p;
    medium.T   = T;
    medium.Xi = X[1:Medium.nXi];
  end FixedAmbient_pTX;
  
  model FixedAmbient_phX 
    "Ambient pressure, specific enthalpy and mass fraction source" 
    extends Interfaces.PartialSource;
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p=
        Medium.reference_p "Ambient pressure";
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy h=
        1.e4 "Ambient specific enthalpy";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
      "Ambient mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nXi > 0));
    annotation (
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
Model <b>FixedAmbient_ph</b> defines constant values for ambient conditions:
</p>
<ul>
<li> Ambient pressure.</li>
<li> Ambient specific enthalpy.</li>
<li> Ambient mass fractions (only for multi-substance flow).</li>
</ul>
<p>
Note, that ambient temperature, specific enthalpy
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
  equation 
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X);
    medium.p = p;
    medium.h = h;
    medium.Xi = X[1:Medium.nXi];
  end FixedAmbient_phX;
  
  annotation (Documentation(info="<html>
<p>
Package <b>Sources</b> contains generic sources for fluid connectors
to define fixed or prescribed ambient conditions.
</p>
</html>"));
  model PrescribedAmbient_pTX 
    "Ambient with prescribed pressure, temperature and composition" 
    extends Interfaces.PartialSource;
    parameter SI.Pressure p = 101325 "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(p_in)==0)));
    parameter SI.Temperature T = Modelica.SIunits.Conversions.from_degC(20) 
      "Fixed value of temperature" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(T_in)==0)));
    parameter SI.MassFraction X[Medium.nX] = Medium.reference_X 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(X_in)==0) or Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in(redeclare type SignalType = 
          SI.Pressure) "Prescribed ambient pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput T_in(
      redeclare type SignalType = SI.Temperature) 
      "Prescribed ambient temperature" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](redeclare type 
        SignalType = SI.MassFraction) "Prescribed ambient composition" 
      annotation (extent=[-140,-80; -100,-40]);
    annotation (
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
Model <b>FixedAmbient_pt</b> defines constant values for ambient conditions:
</p>
<ul>
<li> Prescribed ambient pressure via input signal p.</li>
<li> Prescribed ambient temperature via input signal T.</li>
<li> Fixed ambient mass fractions (only for multi-substance flow).</li>
</ul>
<p>
Note, that ambient temperature
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"),
      Diagram);
  equation 
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      true, X_in);
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
  end PrescribedAmbient_pTX;
  
  model PrescribedAmbient_phX 
    "Ambient with prescribed pressure, specific enthalpy and composition" 
    extends Interfaces.PartialSource;
    parameter SI.Pressure p = 101325 "Fixed value of pressure" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(p_in)==0)));
    parameter SI.SpecificEnthalpy h = 1e4 "Fixed value of specific enthalpy" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(h_in)==0)));
    parameter SI.MassFraction X[Medium.nX] = Medium.reference_X 
      "Fixed value of composition" 
      annotation (Evaluate = true,
                  Dialog(enable = (cardinality(X_in)==0) or Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput p_in(
      redeclare type SignalType = SI.Pressure) "Prescribed ambient pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput h_in(
      redeclare type SignalType = SI.SpecificEnthalpy) 
      "Prescribed ambient specific enthalpy" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = SI.MassFraction) 
      "Prescribed ambient composition" 
      annotation (extent=[-140,-80; -100,-40]);
    annotation (
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
Model <b>PrescribedAmbient_ph</b> defines values for ambient conditions:
</p>
<ul>
<li> Prescribed ambient pressure via input signal p_ambient.</li>
<li> Prescribed ambient specific enthalpy via input signal h_ambient.</li>
<li> Fixed ambient mass fractions (only for multi-substance flow).</li>
</ul>
<p>
Note, that ambient specific enthalpy
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
  equation 
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      true, X_in);
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
  end PrescribedAmbient_phX;
  
  model PrescribedMassFlowRate_TX 
    "Ideal pump that produces a prescribed mass flow with prescribed temperature and mass fraction" 
    extends Interfaces.PartialSource;
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate going out of the fluid port";
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T=
        Modelica.SIunits.Conversions.from_degC(20) 
      "Fixed value of the fluid temperature";
    parameter Medium.MassFraction X[Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
      "Fixed value of the fluid composition" 
      annotation (Dialog(enable = Medium.nXi > 0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in(
      redeclare type SignalType = SI.MassFlowRate) "Prescribed mass flow rate" 
      annotation (extent=[-128,40; -88,80]);
    Modelica.Blocks.Interfaces.RealInput T_in(
      redeclare type SignalType = SI.Temperature) 
      "Prescribed fluid temperature" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = SI.MassFraction) 
      "Prescribed fluid composition" 
      annotation (extent=[-130,-80; -90,-40]);
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(
        Line(points=[-90,-60; -72,-60; -74,-60], style(color=3, rgbcolor={0,0,255})),
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
        Text(extent=[-142, 142; 156, 88], string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Line(points=[-88,60; -74,60], style(color=3, rgbcolor={0,0,255})),
        Text(
          extent=[-148,-90; 156,-134],
          style(color=0),
          string="%m_flow"),
        Text(
          extent=[-158,122; -52,74],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="m_flow"),
        Text(
          extent=[-158,40; -56,0],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="T"),
        Text(
          extent=[-156,-18; -54,-58],
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
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X);
    if cardinality(m_flow_in)==0 then
      m_flow_in = m_flow;
    end if;
    if cardinality(T_in)==0 then
      T_in = T;
    end if;
    if cardinality(X_in)==0 then
      X_in = X;
    end if;
    medium.T = T;
    medium.Xi = X[1:Medium.nXi];
    port.m_flow = -m_flow;
  end PrescribedMassFlowRate_TX;
  
  model PrescribedMassFlowRate_hX 
    "Ideal pump that produces a prescribed mass flow with prescribed specific enthalpy and mass fraction" 
    extends Interfaces.PartialSource;
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate going out of the fluid port";
    parameter Medium.SpecificEnthalpy h = 1e4 
      "Fixed value of the fluid specific enthalpy";
    parameter Medium.MassFraction X[Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
      "Fixed value of the fluid composition" 
      annotation (Dialog(enable=Medium.nXi>0));
    Modelica.Blocks.Interfaces.RealInput m_flow_in(
      redeclare type SignalType = SI.MassFlowRate) "Prescribed mass flow rate" 
      annotation (extent=[-128,40; -88,80]);
    Modelica.Blocks.Interfaces.RealInput h_in(
      redeclare type SignalType = SI.SpecificEnthalpy) 
      "Prescribed fluid specific enthalpy" 
      annotation (extent=[-140,-20; -100,20]);
    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](
      redeclare type SignalType = SI.MassFraction) 
      "Prescribed fluid composition" 
      annotation (extent=[-130,-80; -90,-40]);
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(
        Line(points=[-90,-60; -72,-60; -74,-60], style(color=3, rgbcolor={0,0,255})),
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
        Text(extent=[-142, 142; 156, 88], string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Line(points=[-88,60; -74,60], style(color=3, rgbcolor={0,0,255})),
        Text(
          extent=[-148,-90; 156,-134],
          style(color=0),
          string="%m_flow"),
        Text(
          extent=[-158,122; -52,74],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="m_flow"),
        Text(
          extent=[-158,40; -56,0],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="T"),
        Text(
          extent=[-156,-18; -54,-58],
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
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X);
    if cardinality(m_flow_in)==0 then
      m_flow_in = m_flow;
      port.m_flow = -m_flow;
    else
      port.m_flow = -m_flow_in;
    end if;
    if cardinality(h_in)==0 then
      h_in = h;
    end if;
    if cardinality(X_in)==0 then
      X_in = X;
    end if;
    medium.h = h;
    medium.Xi = X[1:Medium.nXi];
  end PrescribedMassFlowRate_hX;
  
end Sources;
