package Sources "Generic fluid sources" 
  extends Modelica.Icons.Library;
  model FixedAmbient "Ambient source component" 
    extends Interfaces.PartialSource;
    
    parameter Boolean use_p_ambient=true "select p_ambient or d_ambient" 
      annotation (Evaluate=true, Dialog(group=
            "Ambient pressure or ambient density"));
    parameter Modelica_Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=
        101325 "Ambient pressure"          annotation (
       Dialog(group="Ambient pressure or ambient density", enable=use_p_ambient));
    parameter Modelica_Media.Interfaces.PartialMedium.Density d_ambient=1 
      "Ambient density"  annotation (Dialog(group=
            "Ambient pressure or ambient density", enable=not use_p_ambient));
    parameter Boolean use_T_ambient=true "select T_ambient or h_ambient" 
      annotation (Evaluate=true, Dialog(group=
            "Ambient temperature or ambient specific enthalpy"));
    parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) "Ambient temperature" 
      annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                                                                enable=
            use_T_ambient));
    parameter Modelica_Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy" 
      annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                                                                enable=not 
            use_T_ambient));
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nX > 0));
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100, 80; 100, -80], style(
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.incompressible,
      use_p_ambient, Medium.reducedX, Medium.nX, X_ambient);
    if use_p_ambient or Medium.incompressible then
      medium.p = p_ambient;
    else
      medium.d = d_ambient;
    end if;
    
    if use_T_ambient then
      medium.T = T_ambient;
    else
      medium.h = h_ambient;
    end if;
    
    medium.X = X_ambient;
  end FixedAmbient;
  
  model FixedAmbient_pTX 
    "Ambient pressure, temperature and mass fraction source" 
    extends Interfaces.PartialSource;
    
    parameter Modelica_Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=
        101325 "Ambient pressure";
    parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) "Ambient temperature";
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nX > 0));
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100, 80; 100, -80], style(
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.incompressible,
      true, Medium.reducedX, Medium.nX, X_ambient);
    medium.p = p_ambient;
    medium.T = T_ambient;
    medium.X = X_ambient;
  end FixedAmbient_pTX;
  
  model FixedAmbient_phX 
    "Ambient pressure, specific enthalpy and mass fraction source" 
    extends Interfaces.PartialSource;
    
    parameter Modelica_Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=
        101325 "Ambient pressure";
    parameter Modelica_Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy";
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nX > 0));
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100, 80; 100, -80], style(
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.incompressible,
      true, Medium.reducedX, Medium.nX, X_ambient);
      medium.p = p_ambient;
      medium.h = h_ambient;
      medium.X = X_ambient;
  end FixedAmbient_phX;
  
  annotation (Documentation(info="<html>
<p>
Package <b>Sources</b> contains generic sources for fluid connectors
to define fixed or prescribed ambient conditions.
</p>
</html>"));
  model PrescribedAmbient_pT 
    "Prescribed ambient pressure and temperature source with fixed mass fraction" 
    import SI = Modelica.SIunits;
    extends Interfaces.PartialSource;
    
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nX > 0));
    Modelica.Blocks.Interfaces.RealInput p_ambient(redeclare type SignalType = 
          SI.Pressure) "Prescribed ambient pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput T_ambient(redeclare type SignalType = 
          SI.Temperature) "Prescribed ambient temperature" 
      annotation (extent=[-140,-80; -100,-40]);
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100, 80; 100, -80], style(
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
          string="%unit_p"),
        Text(
          extent=[-152,-86; -50,-126],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="%unit_T")),
      Documentation(info="<html>
<p>
Model <b>FixedAmbient_pt</b> defines constant values for ambient conditions:
</p>
<ul>
<li> Prescribed ambient pressure via input signal p_ambient.</li>
<li> Prescribed ambient temperature via input signal T_ambient.</li>
<li> Fixed ambient mass fractions (only for multi-substance flow).</li>
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.incompressible,
      true, Medium.reducedX, Medium.nX, X_ambient);
    medium.p = p_ambient;
    medium.T = T_ambient;
    medium.X = X_ambient;
  end PrescribedAmbient_pT;
  
  model PrescribedAmbient_ph 
    "Prescribed ambient pressure and specific enthalpy source with fixed mass fraction" 
    import SI = Modelica.SIunits;
    extends Interfaces.PartialSource;
    
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m"  annotation (Dialog(group=
            "Only for multi-substance flow", enable=Medium.nX > 0));
    Modelica.Blocks.Interfaces.RealInput p_ambient(redeclare type SignalType = 
          SI.Pressure) "Prescribed ambient pressure" 
      annotation (extent=[-140,40; -100,80]);
    Modelica.Blocks.Interfaces.RealInput h_ambient(redeclare type SignalType = 
          SI.SpecificEnthalpy) "Prescribed ambient specific enthalpy" 
      annotation (extent=[-140,-80; -100,-40]);
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(Ellipse(extent=[-100, 80; 100, -80], style(
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
          string="%unit_p"),
        Text(
          extent=[-152,-86; -50,-126],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="%unit_h")),
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.incompressible,
      true, Medium.reducedX, Medium.nX, X_ambient);
    medium.p = p_ambient;
    medium.h = h_ambient;
    medium.X = X_ambient;
  end PrescribedAmbient_ph;
  
  model FixedMassFlowRate_TX 
    "Ideal pump that produces a constant mass flow rate from a large reservoir at fixed temperature and mass fraction" 
    
    extends Interfaces.PartialSource;
    
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate from an infinite reservoir to the fluid port";
    parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) 
      "Ambient temperature of reservoir";
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m of reservoir" 
      annotation (Dialog(group="Only for multi-substance flow"));
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
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
        Text(extent=[-142, 142; 156, 88], string="%name"),
        Text(
          extent=[-154,-88; 150,-132],
          style(color=0),
          string="%m_flow"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1))),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.incompressible, true,
      Medium.reducedX, Medium.nX, X_ambient);
      medium.T = T_ambient;
      medium.X = X_ambient;
      port.m_flow = -m_flow;
  end FixedMassFlowRate_TX;
  
  model FixedMassFlowRate_hX 
    "Ideal pump that produces a constant mass flow rate from a large reservoir at fixed specific enthalpy and mass fraction" 
    
    extends Interfaces.PartialSource;
    
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate from an infinite reservoir to the fluid port";
    parameter Modelica_Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy  of reservoir";
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m of reservoir" 
      annotation (Dialog(group="Only for multi-substance flow"));
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
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
        Text(extent=[-142, 142; 156, 88], string="%name"),
        Text(
          extent=[-154,-88; 150,-132],
          style(color=0),
          string="%m_flow"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1))),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.incompressible, true,
      Medium.reducedX, Medium.nX, X_ambient);
      medium.h = h_ambient;
      medium.X = X_ambient;
      port.m_flow = -m_flow;
  end FixedMassFlowRate_hX;
  
  model PrescribedMassFlowRate_TX 
    "Ideal pump that produces a prescribed mass flow rate from a large reservoir at fixed temperature and mass fraction" 
    
    extends Interfaces.PartialSource;
    parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) 
      "Ambient temperature of reservoir";
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m of reservoir" 
      annotation (Dialog(group="Only for multi-substance flow"));
    Modelica.Blocks.Interfaces.RealInput m_flow_ambient(redeclare type 
        SignalType =                                                                
          SI.MassFlowRate) "Mass flow rate from large reservoir to fluid port" 
      annotation (extent=[-140,-20; -100,20]);
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
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
        Text(extent=[-154,146; 144,92],   string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Text(
          extent=[-188,-42; -86,-82],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="%unit_m_flow")),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.incompressible, true,
      Medium.reducedX, Medium.nX, X_ambient);
      medium.T = T_ambient;
      medium.X = X_ambient;
    
    port.m_flow = -m_flow_ambient;
  end PrescribedMassFlowRate_TX;
  
  model PrescribedMassFlowRate_hX 
    "Ideal pump that produces a prescribed mass flow rate from a large reservoir at fixed specific enthalpy and mass fraction" 
    
    extends Interfaces.PartialSource;
    parameter Modelica_Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy of reservoir";
    parameter Modelica_Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = zeros(Medium.nX) 
      "Ambient mass fractions m_i/m of reservoir" 
      annotation (Dialog(group="Only for multi-substance flow"));
    Modelica.Blocks.Interfaces.RealInput m_flow_ambient(redeclare type 
        SignalType =                                                                
          SI.MassFlowRate) "Mass flow rate from large reservoir to fluid port" 
      annotation (extent=[-140,-20; -100,20]);
    
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
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
        Text(extent=[-154,146; 144,92],   string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Text(
          extent=[-188,-42; -86,-82],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="%unit_m_flow")),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.incompressible, true,
      Medium.reducedX, Medium.nX, X_ambient);
      medium.h = h_ambient;
      medium.X = X_ambient;
    
    port.m_flow = -m_flow_ambient;
  end PrescribedMassFlowRate_hX;
end Sources;
