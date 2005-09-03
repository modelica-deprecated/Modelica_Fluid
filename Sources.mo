package Sources "Generic fluid sources" 
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  model FixedAmbient "Ambient source component" 
    extends Interfaces.PartialSource;
    
    parameter Boolean use_p_ambient=true "select p_ambient or d_ambient" 
      annotation (Evaluate=true, Dialog(group=
            "Ambient pressure or ambient density"));
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=
        Medium.reference_p "Ambient pressure" 
                                           annotation (
       Dialog(group="Ambient pressure or ambient density", enable=use_p_ambient));
    parameter Modelica.Media.Interfaces.PartialMedium.Density d_ambient=1 
      "Ambient density"  annotation (Dialog(group=
            "Ambient pressure or ambient density", enable=not use_p_ambient));
    parameter Boolean use_T_ambient=true "select T_ambient or h_ambient" 
      annotation (Evaluate=true, Dialog(group=
            "Ambient temperature or ambient specific enthalpy"));
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) "Ambient temperature" 
      annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                                                                enable=
            use_T_ambient));
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy" 
      annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                                                                enable=not 
            use_T_ambient));
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[Medium.nX](
         quantity=Medium.substanceNames)=Medium.reference_X 
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      use_p_ambient, X_ambient);
    
    if use_p_ambient or Medium.singleState then
      medium.p = p_ambient;
    else
      medium.d = d_ambient;
    end if;
    
    if use_T_ambient then
      medium.T = T_ambient;
    else
      medium.h = h_ambient;
    end if;
    
    medium.Xi = X_ambient[1:Medium.nXi];
  end FixedAmbient;
  
  model FixedAmbient_pTX 
    "Ambient pressure, temperature and mass fraction source" 
    extends Interfaces.PartialSource;
    
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=
        Medium.reference_p "Ambient pressure";
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) "Ambient temperature";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[Medium.nX](
         quantity=Medium.substanceNames) = Medium.reference_X 
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      true, X_ambient);
    medium.p   = p_ambient;
    medium.T   = T_ambient;
    medium.Xi = X_ambient[1:Medium.nXi];
  end FixedAmbient_pTX;
  
  model FixedAmbient_phX 
    "Ambient pressure, specific enthalpy and mass fraction source" 
    extends Interfaces.PartialSource;
    
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure 
      p_ambient=
        Medium.reference_p "Ambient pressure";
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
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
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X_ambient);
      medium.p = p_ambient;
      medium.h = h_ambient;
      medium.Xi = X_ambient[1:Medium.nXi];
  end FixedAmbient_phX;
  
  annotation (Documentation(info="<html>
<p>
Package <b>Sources</b> contains generic sources for fluid connectors
to define fixed or prescribed ambient conditions.
</p>
</html>"));
  model PrescribedAmbient_pT 
    "Prescribed ambient pressure and temperature source with fixed mass fraction" 
    extends Interfaces.PartialSource;
    
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
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
          string="p"),
        Text(
          extent=[-152,-86; -50,-126],
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
</html>"),
      Diagram);
    
  equation 
    Modelica_Fluid.Utilities.checkAmbient(Medium.mediumName, Medium.singleState,
      true, X_ambient);
    medium.p = p_ambient;
    medium.T = T_ambient;
    medium.Xi = X_ambient[1:Medium.nXi];
  end PrescribedAmbient_pT;
  
  model PrescribedAmbient_ph 
    "Prescribed ambient pressure and specific enthalpy source with fixed mass fraction" 
    extends Interfaces.PartialSource;
    
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
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
      true, X_ambient);
    medium.p = p_ambient;
    medium.h = h_ambient;
    medium.Xi = X_ambient[1:Medium.nXi];
  end PrescribedAmbient_ph;
  
  model FixedMassFlowRate_TX 
    "Ideal pump that produces a constant mass flow rate from a large reservoir at fixed temperature and mass fraction" 
    
    extends Interfaces.PartialSource;
    
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate from an infinite reservoir to the fluid port";
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) 
      "Ambient temperature of reservoir";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
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
    Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X_ambient);
      medium.T = T_ambient;
      medium.Xi = X_ambient[1:Medium.nXi];
      port.m_flow = -m_flow;
  end FixedMassFlowRate_TX;
  
  model FixedMassFlowRate_hX 
    "Ideal pump that produces a constant mass flow rate from a large reservoir at fixed specific enthalpy and mass fraction" 
    
    extends Interfaces.PartialSource;
    
    parameter Medium.MassFlowRate m_flow 
      "Fixed mass flow rate from an infinite reservoir to the fluid port";
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy  of reservoir";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
      "Ambient mass fractions m_i/m of reservoir, reduced" 
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
    Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X_ambient);
      medium.h = h_ambient;
      medium.Xi = X_ambient[1:Medium.nXi];
      port.m_flow = -m_flow;
  end FixedMassFlowRate_hX;
  
  model PrescribedMassFlowRate_TX 
    "Ideal pump that produces a prescribed mass flow rate from a large reservoir at fixed temperature and mass fraction" 
    
    extends Interfaces.PartialSource;
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=
        Modelica.SIunits.Conversions.from_degC(20) 
      "Ambient temperature of reservoir";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
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
          string="m_flow")),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X_ambient);
    medium.T = T_ambient;
    medium.Xi = X_ambient[1:Medium.nXi];
    port.m_flow = -m_flow_ambient;
  end PrescribedMassFlowRate_TX;
  
  model PrescribedMassFlowRate_hX 
    "Ideal pump that produces a prescribed mass flow rate from a large reservoir at fixed specific enthalpy and mass fraction" 
    
    extends Interfaces.PartialSource;
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy 
      h_ambient=
        1.e4 "Ambient specific enthalpy of reservoir";
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_ambient[
      Medium.nX](quantity=Medium.substanceNames) = Medium.reference_X 
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
        Text(extent=[-150,138; 148,84],   string="%name"),
        Ellipse(extent=[-26, 30; -18, 22], style(color=1, fillColor=1)),
        Text(
          extent=[-188,-42; -86,-82],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255}),
          string="m_flow")),
      Window(
        x=0.45,
        y=0.01,
        width=0.44,
        height=0.65),
      Diagram);
  equation 
    Utilities.checkAmbient(Medium.mediumName, Medium.singleState, true, X_ambient);
    
    medium.h = h_ambient;
    medium.Xi = X_ambient[1:Medium.nXi];
    port.m_flow = -m_flow_ambient;
  end PrescribedMassFlowRate_hX;
  
  model SourceP "Prescribed pressure ambient" 
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
    Medium.BaseProperties medium(p(start=p0),T(start=T),Xi(start=Xnom[1:Medium.nXi]));
    parameter SI.Pressure p0=101325 "Nominal pressure";
    parameter Real R=0 "Hydraulic resistance";
    parameter Medium.Temperature T=300 "Nominal temperature";
    parameter Medium.MassFraction Xnom[Medium.nX]=Medium.reference_X 
      "Nominal medium composition";
    
    Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[80,-20; 120,20]);
    Modelica.Blocks.Interfaces.RealInput in_p 
      annotation (extent=[-70,54; -50,74], rotation=-90);
    Modelica.Blocks.Interfaces.RealInput in_T 
      annotation (extent=[-10,80; 10,100], rotation=270);
    Modelica.Blocks.Interfaces.RealInput in_X[Medium.nX] 
      annotation (extent=[50,52; 70,72], rotation=270);
  equation 
    if R == 0 then
      port.p = medium.p;
    else
      port.p = medium.p + port.m_flow*R;
    end if;
    
    medium.p = in_p;
    if cardinality(in_p)==0 then
      in_p = p0 "Pressure set by parameter";
    end if;
    
    medium.T = in_T;
    if cardinality(in_T)==0 then
      in_T = T "Temperature set by parameter";
    end if;
    
    medium.Xi = in_X[1:Medium.nXi];
    if cardinality(in_X)==0 then
      in_X = Xnom "Composition set by parameter";
    end if;
    
    port.H_flow = semiLinear(port.m_flow,port.h,medium.h);
    port.mXi_flow = semiLinear(port.m_flow,port.Xi,medium.Xi);
    
    annotation (Icon(
           Ellipse(extent=[-80,80; 80,-80],     style(
            color=69,
            gradient=3,
            fillColor=69))),
                      Diagram,
      Documentation(info="<html>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package.In the case of multiple componet, variable composition gases, the nominal gas composition is given by <tt>Xnom</tt>, whose default value is <tt>Medium.reference_X</tt> .
<p>If <tt>R</tt> is set to zero, the pressure source is ideal; otherwise, the outlet pressure decreases proportionally to the outgoing flowrate.</p>
<p>If the <tt>in_p</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>p0</tt>.</p>
<p>If the <tt>in_T</tt> connector is wired, then the source temperature is given by the corresponding signal, otherwise it is fixed to <tt>T</tt>.</p>
<p>If the <tt>in_X</tt> connector is wired, then the source massfraction is given by the corresponding signal, otherwise it is fixed to <tt>Xnom</tt>.</p>
</html>", revisions="<html>
<ul>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>p0fix</tt> and <tt>Tfix</tt> and <tt>Xfix</tt>; the connection of external signals is now detected automatically.</li> <br> Adapted to Modelica.Media
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>
"));
  end SourceP;
  
  model SourceFlow "Flowrate source for medium flows" 
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
    Medium.BaseProperties medium(p(start=p0),T(start=T),Xi(start=Xnom[1:Medium.nXi]));
    parameter SI.Pressure p0=101325 "Nominal pressure";
    parameter SI.Temperature T=300 "Nominal temperature";
    parameter SI.MassFraction Xnom[Medium.nX]=Medium.reference_X 
      "Nominal medium composition";
    parameter SI.MassFlowRate m_flow0=0 "Nominal mass flowrate";
    parameter Real.HydraulicConductance G(unit="kg/(s.Pa)")=0 
      "HydraulicConductance";
    
    SI.MassFlowRate m_flow;
    
    Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[80,-20; 120,20]);
    Modelica.Blocks.Interfaces.RealInput in_m_flow0 
      annotation (extent=[-86,48; -60,74], rotation=-90);
    Modelica.Blocks.Interfaces.RealInput in_T 
      annotation (extent=[-6,98; -32,68], rotation=-270);
    Modelica.Blocks.Interfaces.RealInput in_X[Medium.nX] 
      annotation (extent=[14,54; 40,86], rotation=-90);
  equation 
    
    if G == 0 then
      port.m_flow = -m_flow;
    else
      port.m_flow = -m_flow + (port.p - p0)*G;
    end if;
    
    m_flow0 = in_m_flow0;
    if cardinality(in_m_flow0)==0 then
      in_m_flow0 = m_flow0 "Flow rate set by parameter";
    end if;
    
    medium.T = in_T;
    if cardinality(in_T)==0 then
      in_T = T "Temperature set by parameter";
    end if;
    
    medium.Xi = in_X[1:Medium.nXi];
    if cardinality(in_X)==0 then
      in_X = Xnom "Composition set by parameter";
    end if;
    
    port.p = medium.p;
    port.H_flow = semiLinear(port.m_flow,port.h,medium.h);
    port.mXi_flow = semiLinear(port.m_flow,port.Xi,medium.Xi);
    
    annotation (Icon(
        Rectangle(extent=[20,48; 100,-72],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[38,28; 100,-52],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[-100,68; 60,-92],   style(fillColor=7)),
        Polygon(points=[-60,58; 60,-12; -60,-80; -60,58],   style(color=73,
              fillColor=73)),
        Text(
          extent=[-54,20; 16,-42],
          style(color=41, fillColor=41),
          string="m"),
        Text(extent=[-166,-88; 132,-142], string="%name"),
        Ellipse(extent=[-26,18; -18,10],   style(color=1, fillColor=1))),
                      uses(Modelica(version="1.6")),
      Documentation(info="<html>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the nominal gas composition is given by <tt>Xnom</tt>,whose default value is <tt>Medium.reference_X</tt> .
<p>If <tt>G</tt> is set to zero, the flowrate source is ideal; otherwise, the outgoing flowrate decreases proportionally to the outlet pressure.</p>
<p>If the <tt>in_w0</tt> connector is wired, then the source massflowrate is given by the corresponding signal, otherwise it is fixed to <tt>w0</tt>.</p>
<p>If the <tt>in_T</tt> connector is wired, then the source temperature is given by the corresponding signal, otherwise it is fixed to <tt>T</tt>.</p>
<p>If the <tt>in_X</tt> connector is wired, then the source massfraction is given by the corresponding signal, otherwise it is fixed to <tt>Xnom</tt>.</p>
</html>", revisions="<html>
<ul>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>w0fix</tt> and <tt>Tfix</tt> and <tt>Xfix</tt>; the connection of external signals is now detected automatically.</li> <br> Adapted to Modelica.Media
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"),
      Diagram);
  end SourceFlow;
end Sources;
