package PressureLosses 
  "Models and functions providing pressure loss correlations " 
model StaticHead 
    "Models the static head between two ports at different heights" 
  extends BaseClasses.Common.PartialTwoPortTransport;
  parameter SI.Length height_ab "Height(port_b) - Height(port_a)";
  Medium.Density d "Fluid density";
  outer Components.Ambient ambient "Ambient conditions";
  annotation (Icon(
      Rectangle(extent=[-100,60; 100,-60],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,34; 100,-36],   style(
          color=69,
          gradient=2,
          fillColor=69)),
      Text(
        extent=[-150,140; 150,80],
        string="%name",
        style(gradient=2, fillColor=69))), Documentation(info="<html>
<p>
This model describes the static head due to the relative height between the two connectors. No mass, energy and momentum storage, and no pressure drop due to friction are considered.
</p>
</html>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added to Modelica_Fluid</li>
</ul>
</html>"));
equation 
 d = if dp > 0 then medium_a.d else medium_b.d;
 dp = height_ab*ambient.g*d;
end StaticHead;
  
model SimpleGenericOrifice 
    "Simple generic orifice defined by pressure loss coefficient and diameter (only for flow from port_a to port_b)" 
    import SI = Modelica.SIunits;
    extends BaseClasses.Common.PartialGuessValueParameters;
    extends BaseClasses.Common.PartialTwoPortTransport(
                medium_a(p(start=p_start), h(start=h_start),
                         T(start=T_start), Xi(start=X_start[1:Medium.nXi])),
                medium_b(p(start=p_start), h(start=h_start),
                         T(start=T_start), Xi(start=X_start[1:Medium.nXi])));
    
  parameter Real zeta "Loss factor for flow of port_a -> port_b";
  parameter SI.Diameter diameter 
      "Diameter at which zeta is defined (either port_a or port_b)";
  parameter Boolean from_dp = true 
      "= true, use m_flow = f(dp) else dp = f(m_flow)" 
    annotation (Evaluate=true, Dialog(tab="Advanced"));
  parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow if |dp| >= dp_small" 
    annotation(Dialog(tab="Advanced", enable=from_dp));
  parameter SI.MassFlowRate m_flow_small = 0.01 
      "Turbulent flow if |m_flow| >= m_flow_small" 
    annotation(Dialog(tab="Advanced", enable=not from_dp));
  annotation (
    preferedView="info",
    Diagram,
    Icon(
      Text(
        extent=[-148,109; 148,58],
        string="%name",
        style(gradient=2, fillColor=69)),
      Line(points=[-60, -50; -60, 50; 60, -50; 60, 50], style(color=0,
            thickness=2)),
      Line(points=[-60, 0; -100, 0], style(color=69)),
      Line(points=[60, 0; 100, 0], style(color=69)),
        Text(
          extent=[-168,-92; 180,-134],
          string="zeta=%zeta",
          style(color=0, rgbcolor={0,0,0})),
        Line(points=[-50,-70; 50,-70], style(color=69, rgbcolor={0,128,255})),
        Polygon(points=[24,-60; 24,-80; 50,-70; 24,-60], style(
            color=69,
            rgbcolor={0,128,255},
            fillColor=69,
            rgbfillColor={0,128,255}))),
    Documentation(info="<html>
<p>
This pressure drop component defines a
simple, generic orifice, where the loss factor &zeta; is provided
for one flow direction (e.g., from loss table of a book):
</p>
 
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 8*&zeta;/(&pi;^2*D^4*&rho;) * m_flow*|m_flow|
</pre>
 
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> D is the diameter of the orifice at the position where
     &zeta; is defined (either at port_a or port_b). If the orifice has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> &zeta; is the loss factor with respect to D 
     that depends on the geometry of
     the orifice. In the turbulent flow regime, it is assumed that
     &zeta; is constant.<br>
     For small mass flow rates, the flow is laminar and is approximated 
     by a polynomial that has a finite derivative for m_flow=0.</li>
<li> v is the mean velocity.</li>
<li> &rho; is the upstream density.</li>
</ul>
 
<p>
Since the pressure loss factor zeta is provided only for a mass flow
from port_a to port_b, the pressure loss is not correct when the
flow is reversing. If reversing flow only occurs in a short time interval,
this is most likely uncritical. If significant reversing flow
can appear, this component should not be used.
</p>
</html>"),
      Coordsys(grid=[1,1], scale=0));
equation 
  if from_dp then
     m_flow = BaseClasses.PressureLosses.SimpleGenericOrifice.massFlowRate_dp(dp, medium_a.d, medium_b.d, diameter, zeta);
  else
     dp = BaseClasses.PressureLosses.SimpleGenericOrifice.pressureLoss_m_flow(m_flow, medium_a.d, medium_b.d, diameter, zeta);
  end if;
end SimpleGenericOrifice;
  
  model WallFrictionAndGravity 
    "Pressure drop in pipe due to wall friction and gravity (for both flow directions)" 
    import SI = Modelica.SIunits;
    
    extends BaseClasses.Common.PartialGuessValueParameters;
    extends BaseClasses.Common.PartialTwoPortTransport(
                  medium_a(p(start=p_start), h(start=h_start),
                           T(start=T_start), Xi(start=X_start[1:Medium.nXi])),
                  medium_b(p(start=p_start), h(start=h_start),
                           T(start=T_start), Xi(start=X_start[1:Medium.nXi])));
    
    replaceable package WallFriction = 
      BaseClasses.PressureLosses.WallFriction.QuadraticTurbulent 
      extends BaseClasses.PressureLosses.WallFriction.PartialWallFriction 
      "Characteristic of wall friction"  annotation(choicesAllMatching=true);
    
    parameter SI.Length length "Length of pipe";
    parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    parameter SI.Length height_ab = 0.0 "Height(port_b) - Height(port_a)" 
                                                                       annotation(Evaluate=true);
    parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
        annotation(Dialog(enable=WallFriction.use_roughness));
    
    parameter Boolean use_nominal = false 
      "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
                                                                                                    annotation(Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal = 0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                              annotation(Dialog(enable=use_nominal));
    parameter SI.Density d_nominal = 0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(enable=use_nominal));
    
    parameter Boolean show_Re = false 
      "= true, if Reynolds number is included for plotting" 
       annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter Boolean from_dp=true 
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
      annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
      annotation(Dialog(tab="Advanced", enable=from_dp and WallFriction.use_dp_small));
    parameter SI.MassFlowRate m_flow_small = 0.01 
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
      annotation(Dialog(tab="Advanced", enable=not from_dp and WallFriction.use_m_flow_small));
    SI.ReynoldsNumber Re = Utilities.ReynoldsNumber_m_flow(m_flow, (eta_a+eta_b)/2, diameter) if show_Re 
      "Reynolds number of pipe";
    
    outer Components.Ambient ambient "Ambient conditions";
    
    annotation (defaultComponentName="pipe",Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-152,119; 138,70],
          string="%name",
          style(gradient=2, fillColor=69))), Documentation(info="<html>
<p>
This model describes pressure losses due to <b>wall friction</b> in a pipe
and due to gravity.
It is assumed that no mass or energy is stored in the pipe. 
Correlations of different complexity and validity can be
seleted via the replaceable package <b>WallFriction</b> (see parameter menu below).
The details of the pipe wall friction model are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
Basically, different variants of the equation
</p>
 
<pre>
   dp = &lambda;(Re,<font face=\"Symbol\">D</font>)*(L/D)*&rho;*v*|v|/2
</pre>
 
<p>
are used, where the friction loss factor &lambda; is shown
in the next figure:
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
<p>
By default, the correlations are computed with media data
at the actual time instant.
In order to reduce non-linear equation systems, parameter
<b>use_nominal</b> provides the option
to compute the correlations with constant media values
at the desired operating point. This might speed-up the
simulation and/or might give a more robust simulation.
</p>
</html>"),
      Diagram(
        Rectangle(extent=[-100,64; 100,-64], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=8)),
        Rectangle(extent=[-100,50; 100,-49], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
        Line(points=[-60,-49; -60,50], style(
            color=3,
            rgbcolor={0,0,255},
            arrow=3,
            fillColor=3,
            rgbfillColor={0,0,255},
            fillPattern=1)),
        Text(
          extent=[-50,16; 6,-10],
          style(
            color=3,
            rgbcolor={0,0,255},
            arrow=3,
            fillColor=3,
            rgbfillColor={0,0,255},
            fillPattern=1),
          string="diameter"),
        Line(points=[-100,74; 100,74], style(
            color=3,
            rgbcolor={0,0,255},
            arrow=3,
            fillColor=3,
            rgbfillColor={0,0,255},
            fillPattern=1)),
        Text(
          extent=[-34,92; 34,74],
          style(
            color=3,
            rgbcolor={0,0,255},
            arrow=3,
            fillColor=3,
            rgbfillColor={0,0,255},
            fillPattern=1),
          string="length")),
      Coordsys(grid=[1,1], scale=0));
  protected 
    SI.DynamicViscosity eta_a = if not WallFriction.use_eta then 1.e-10 else 
                                (if use_nominal then eta_nominal else Medium.dynamicViscosity(medium_a));
    SI.DynamicViscosity eta_b = if not WallFriction.use_eta then 1.e-10 else 
                                (if use_nominal then eta_nominal else Medium.dynamicViscosity(medium_b));
    SI.Density d_a = if use_nominal then d_nominal else medium_a.d;
    SI.Density d_b = if use_nominal then d_nominal else medium_b.d;
  equation 
    if from_dp and not WallFriction.dp_is_zero then
       m_flow = WallFriction.massFlowRate_dp(dp-height_ab*ambient.g*(d_a+d_b)/2,
                                             d_a, d_b, eta_a, eta_b, length, diameter, roughness, dp_small);
    else
       dp = WallFriction.pressureLoss_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, roughness, m_flow_small)
            + height_ab*ambient.g*(d_a+d_b)/2;
    end if;
  end WallFrictionAndGravity;
  extends Modelica.Icons.Library;
  
model SuddenExpansion 
    "Pressure drop in pipe due to suddenly expanding area (for both flow directions)" 
    import SI = Modelica.SIunits;
  extends BaseClasses.PressureLosses.QuadraticTurbulent.BaseModel(
     final data = BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.suddenExpansion(D_a, D_b));
  parameter SI.Diameter D_a "Inner diameter of pipe at port_a";
  parameter SI.Diameter D_b "Inner diameter of pipe at port_b";
    
  annotation (
    Diagram(
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
        string="D_a"),
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
        string="D_b")),
    Icon(
      Text(
        extent=[-132,144; 138,98],
        string="%name",
        style(gradient=2, fillColor=69)),
      Rectangle(extent=[-100,80; 100,-80],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,20; 0,-20],     style(
          color=69,
          gradient=2,
          fillColor=69)),
      Rectangle(extent=[0,60; 100,-60],      style(
          color=69,
          gradient=2,
          fillColor=69))),
      Documentation(info="<html>
 
</html>"));
end SuddenExpansion;
  
model SharpEdgedOrifice 
    "Pressure drop due to sharp edged orifice (for both flow directions)" 
    import SI = Modelica.SIunits;
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
  extends BaseClasses.PressureLosses.QuadraticTurbulent.BaseModel(
     final data = BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.sharpEdgedOrifice(D_pipe, D_min, L, alpha));
  parameter SI.Diameter D_pipe 
      "Inner diameter of pipe (= same at port_a and port_b)";
  parameter SI.Diameter D_min "Smallest diameter of orifice";
  parameter SI.Diameter L "Length of orifice";
  parameter NonSI.Angle_deg alpha "Angle of orifice";
  annotation (defaultComponentName="orifice",
    Documentation(info="<html>
</html>"),
    Icon(
      Text(
        extent=[-148,136; 148,92],
        string="%name",
        style(gradient=2, fillColor=69)),
      Rectangle(extent=[-100,80; 100,-80],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,60; 100,-60],   style(
          color=69,
          gradient=2,
          fillColor=69)),
      Polygon(points=[-24,60; -24,12; 36,50; 36,60; -24,60], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8)),
      Polygon(points=[-22,-10; -22,-60; 38,-60; 38,-50; -22,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8))),
    Diagram(       Rectangle(extent=[-100,60; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Polygon(points=[-30,60; -30,12; 30,50; 30,60; -30,60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8)),
      Polygon(points=[-30,-10; -30,-60; 30,-60; 30,-50; -30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8)),
      Line(points=[-82,-60; -82,60], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-78,16; -44,-8],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="D_pipe"),
      Line(points=[-30,-10; -30,12], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-24,14; 8,-10],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="D_min"),
      Text(
        extent=[-20,84; 18,70],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="L"),
      Line(points=[30,68; -30,68],   style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Line(points=[16,40; 32,18; 36,-2; 34,-20; 20,-42], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8)),
      Text(
        extent=[38,8; 92,-6],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8),
        string="alpha")));
end SharpEdgedOrifice;
  
  annotation (Documentation(info="<html>
<p>
This sublibrary contains models and functions providing pressure 
loss correlations. All models in this library have the property
that no mass and no energy is stored in the component. Therefore,
none of the models has a state. The basic correlations are implemented
with functions of sublibrary
<a href=\"Modelica://Modelica_Fluid.PressureLosses.Utilities\">PressureLosses.Utilities</a>
These functions might also be directly called 
(e.g. in another component implementation).
</p>
 
<p>
All functions are continuous and have a finite, non-zero, smooth, first derivative.
The functions are all guaranteed to be strict monontonically increasing.
The mentioned properties guarantee that a unique inverse of every
function exists. Note, the usual quadratic pressure loss correlation
</p>
 
<ul>
<li> in the form m_flow = f(dp) has an infinite derivative at zero 
     mass flow rate and is therefore problematic to use.</li>
<li> in the form dp = f(m_flow) has a zero derivative at zero mass flow rate
     and is therefore problematic to invert, since the inverse function has
     then an infinite derivative at zero mass flow rate.</li>
</ul>
<p>
The two mentioned problems are solved in this package by approximating
the characteristics around zero mass flow rates with appropriate
polynomials. The monotonicity is guaranteed using results from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
 
</html>", revisions="<html>
<ul>
<li><i>Jan. 3, 2006</i>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    New design and implementation based on previous iterations.</li>
</ul>
</html>"));
  
model PressureDropPipe 
    "Obsolet component (use instead PressureLosses.WallFrictionAndGravity)" 
  extends Modelica_Fluid.BaseClasses.Common.PartialTwoPortTransport;
  extends Modelica_Fluid.BaseClasses.PressureLosses.PipeFriction;
  annotation (Icon(
      Rectangle(extent=[-100,60; 100,-60],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,34; 100,-36],   style(
          color=69,
          gradient=2,
          fillColor=69)),
      Text(
        extent=[-150,140; 150,80],
        string="%name",
        style(gradient=2, fillColor=69))), Documentation(info="<html>
<p>
This model describes pressure losses due to friction in a pipe. It is assumed that no mass or energy is stored in the pipe. 
The details of the pipe friction model are described
<a href=\"Modelica://Modelica_Fluid.Utilities.PipeFriction\">here</a>.
</p>
</html>"));
equation 
  if frictionType == Modelica_Fluid.Types.FrictionTypes.DetailedFriction then
     if from_dp then
        d = if dp >= 0 then medium_a.d else medium_b.d;
        eta = if dp >= 0 then Medium.dynamicViscosity(medium_a) else 
                             Medium.dynamicViscosity(medium_b);
     else
        d = if m_flow >= 0 then medium_a.d else medium_b.d;
        eta = if m_flow >= 0 then Medium.dynamicViscosity(medium_a) else 
                             Medium.dynamicViscosity(medium_b);
     end if;
  else
    // Assign dummy values for auxiliary variables
     d = 0;
     eta = 0;
  end if;
end PressureDropPipe;
end PressureLosses;
