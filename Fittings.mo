within Modelica_Fluid;
package Fittings "Adapt connections of fluid models"
     extends Modelica_Fluid.Icons.VariantLibrary;

  model TJunctionIdeal
    "Splitting/joining component with static balances for an infinitesimal control volume"
    extends BaseClasses.PartialTJunction;

  equation
    connect(port_1, port_2) annotation (Line(
        points={{-100,0},{100,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(port_1, port_3) annotation (Line(
        points={{-100,0},{0,0},{0,100}},
        color={0,127,255},
        smooth=Smooth.None));

    annotation(Documentation(info="<html>
  This model is the simplest implementation for a splitting/joining component for
  three flows. Its use is not required. It just formulates the balance
  equations in the same way that the connect symmantics would formulate them anyways.
  The main advantage of using this component is, that the user does not get
  confused when looking at the specific enthalpy at each port which might be confusing
  when not using a splitting/joining component. The reason for the confusion is that one exmanins the mixing
  enthalpy of the infinitesimal control volume introduced with the connect statement when
  looking at the specific enthalpy in the connector which
  might not be equal to the specific enthalpy at the port in the \"real world\".</html>"),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
  end TJunctionIdeal;

  annotation (Documentation(info="<html>
 
</html>"));
  model TJunctionVolume
    "Splitting/joining component with static balances for a dynamic control volume"
    extends BaseClasses.PartialTJunction;
    extends Modelica_Fluid.Vessels.BaseClasses.PartialLumpedVolume(
                                                    Qs_flow = 0);

    parameter SI.Volume V "Mixing volume inside junction";

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
    assert(cardinality(port_1) <= 1,"
port_1 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
    assert(cardinality(port_2) <= 1,"
port_2 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
    assert(cardinality(port_3) <= 1,"
port_3 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");

    // Boundary conditions
    port_1.h_outflow = medium.h;
    port_2.h_outflow = medium.h;
    port_3.h_outflow = medium.h;

    port_1.Xi_outflow = medium.Xi;
    port_2.Xi_outflow = medium.Xi;
    port_3.Xi_outflow = medium.Xi;

    port_1.C_outflow = C;
    port_2.C_outflow = C;
    port_3.C_outflow = C;

    // Mass balances
    fluidVolume = V;
    ms_flow = port_1.m_flow + port_2.m_flow + port_3.m_flow "Mass balance";
    msXi_flow = port_1.m_flow*actualStream(port_1.Xi_outflow)
                + port_2.m_flow*actualStream(port_2.Xi_outflow)
                + port_3.m_flow*actualStream(port_3.Xi_outflow)
      "Component mass balances";

    msC_flow  = port_1.m_flow*actualStream(port_1.C_outflow)
              + port_2.m_flow*actualStream(port_2.C_outflow)
              + port_3.m_flow*actualStream(port_3.C_outflow)
      "Trace substance mass balances";

    // Momentum balance (suitable for compressible media)
    port_1.p = medium.p;
    port_2.p = medium.p;
    port_3.p = medium.p;

    // Energy balance
    Hs_flow = port_1.m_flow*actualStream(port_1.h_outflow)
              + port_2.m_flow*actualStream(port_2.h_outflow)
              + port_3.m_flow*actualStream(port_3.h_outflow);
    Ws_flow = 0;

    annotation (Documentation(info="<html>
  This model introduces a mixing volume into a junction. 
  This might be useful to examine the non-ideal mixing taking place in a real junction.</html>"),
  Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Ellipse(
            extent={{-9,10},{11,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
  end TJunctionVolume;

  model MultiPort
    "Multiply a port; useful if multiple connections shall be made to a port exposing a state"

    function positiveMax
      input Real x;
      output Real y;
    algorithm
      y :=max(x, 1e-10);
    end positiveMax;

    import Modelica.Constants;

    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);

    // Ports
    parameter Integer nPorts_b=1
      "Number of outlet ports (mass is distributed evenly between the outlet ports";
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
      redeclare package Medium=Medium) 
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.FluidPorts_b[nPorts_b] ports_b(
      redeclare each package Medium=Medium) 
      annotation (Placement(transformation(extent={{30,40},{50,-40}},
                                  rotation=0)));

    Medium.MassFraction[nPorts_b,Medium.nXi] ports_b_Xi_inStream
      "inStream mass fractions at ports_b";
    Medium.ExtraProperty[nPorts_b,Medium.nC] ports_b_C_inStream
      "inStream extra properties at ports_b";

    annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-40,
              -100},{40,100}}), graphics={
          Line(
            points={{-40,0},{40,0}},
            color={0,128,255},
            thickness=1),
          Line(
            points={{-40,0},{40,26}},
            color={0,128,255},
            thickness=1),
          Line(
            points={{-40,0},{40,-26}},
            color={0,128,255},
            thickness=1),
          Text(
            extent={{-150,100},{150,60}},
            lineColor={0,0,255},
            textString="%name")}),
                            Documentation(info="<html>
<p>
This model is useful if multiple connections shall be made to a port of a volume model exposing a state,
like a pipe with ModelStructure av_vb. 
The mixing is shifted into the volume connected to port_a and the result is propageted back to each ports_b.
</p>
<p>
If multiple connections were directly made to the volume,
then ideal mixing would take place in the connection set, outside the volume. This is normally not intended.
</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-40,-100},{
              40,100}}),
              graphics));

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
    for i in 1:nPorts_b loop
      assert(cardinality(ports_b[i]) <= 1,"
each ports_b[i] of boundary shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts_b to add an additional port.
");
    end for;

    // mass and momentum balance
    0 = port_a.m_flow + sum(ports_b.m_flow);
    ports_b.p = fill(port_a.p, nPorts_b);

    // expose stream values from port_a to ports_b
    ports_b.h_outflow = fill(inStream(port_a.h_outflow), nPorts_b);
    ports_b.Xi_outflow = fill(inStream(port_a.Xi_outflow), nPorts_b);
    ports_b.C_outflow = fill(inStream(port_a.C_outflow), nPorts_b);

    // mixing at port_a
    port_a.h_outflow = sum({positiveMax(ports_b[j].m_flow)*inStream(ports_b[j].h_outflow) for j in 1:nPorts_b})
                         / sum({positiveMax(ports_b[j].m_flow) for j in 1:nPorts_b});
    for j in 1:nPorts_b loop
       ports_b_Xi_inStream[j,:] = inStream(ports_b[j].Xi_outflow);
       ports_b_C_inStream[j,:] = inStream(ports_b[j].C_outflow);
    end for;
    for i in 1:Medium.nXi loop
      port_a.Xi_outflow[i] = (positiveMax(ports_b.m_flow)*ports_b_Xi_inStream[:,i])
                           / sum(positiveMax(ports_b.m_flow));
    end for;
    for i in 1:Medium.nC loop
      port_a.C_outflow[i] = (positiveMax(ports_b.m_flow)*ports_b_C_inStream[:,i])
                           / sum(positiveMax(ports_b.m_flow));
    end for;
  end MultiPort;

model SimpleGenericOrifice
    "Simple generic orifice defined by pressure loss coefficient and diameter (only for flow from port_a to port_b)"
      extends Modelica_Fluid.Fittings.BaseClasses.PartialTwoPortTransport;

  parameter Real zeta "Loss factor for flow of port_a -> port_b";
  parameter SI.Diameter diameter
      "Diameter at which zeta is defined (either port_a or port_b)";
  parameter Boolean from_dp = true
      "= true, use m_flow = f(dp) else dp = f(m_flow)" 
    annotation (Evaluate=true, Dialog(tab="Advanced"));
  annotation (
    Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
    Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Text(
            extent={{-150,60},{150,100}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(
            points={{-60,-50},{-60,50},{60,-50},{60,50}},
            color={0,0,0},
            thickness=0.5),
          Line(points={{-60,0},{-100,0}}, color={0,127,255}),
          Line(points={{60,0},{100,0}}, color={0,127,255}),
          Text(
            extent={{-168,-96},{180,-138}},
            lineColor={0,0,0},
            textString="zeta=%zeta")}),
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
</html>"));
equation
  if from_dp then
      m_flow = BaseClasses.SimpleGenericOrifice.massFlowRate_dp(
        dp,
        port_a_d_inflow,
        port_b_d_inflow,
        diameter,
        zeta);
  else
      dp = BaseClasses.SimpleGenericOrifice.pressureLoss_m_flow(
        m_flow,
        port_a_d_inflow,
        port_b_d_inflow,
        diameter,
        zeta);
  end if;
end SimpleGenericOrifice;

model SharpEdgedOrifice
    "Pressure drop due to sharp edged orifice (for both flow directions)"
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
  extends BaseClasses.QuadraticTurbulent.BaseModel(final data=
          BaseClasses.QuadraticTurbulent.LossFactorData.sharpEdgedOrifice(
          diameter,
          leastDiameter,
          length,
          alpha));
  parameter SI.Length length "Length of orifice";
  parameter SI.Diameter diameter
      "Inner diameter of pipe (= same at port_a and port_b)";
  parameter SI.Diameter leastDiameter "Smallest diameter of orifice";
  parameter NonSI.Angle_deg alpha "Angle of orifice";
  annotation (defaultComponentName="orifice",
    Documentation(info="<html>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}), graphics={
          Text(
            extent={{-150,90},{150,130}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Rectangle(
            extent={{-100,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,50},{100,-50}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Polygon(
            points={{-25,50},{-25,7},{35,45},{35,50},{-25,50}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={255,255,255}),
          Polygon(
            points={{-24.5,-5},{-24.5,-50},{35.5,-50},{35.5,-45},{-24.5,-5}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward)}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,60},{100,-60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-30,60},{-30,12},{30,50},{30,60},{-30,60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward),
          Polygon(
            points={{-30,-10},{-30,-60},{30,-60},{30,-50},{-30,-10}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward),
          Line(
            points={{-82,-60},{-82,60}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{-78,16},{-44,-8}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="diameter"),
          Line(
            points={{-30,-10},{-30,12}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{-24,14},{8,-10}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="leastDiameter"),
          Text(
            extent={{-20,84},{18,70}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="length"),
          Line(
            points={{30,68},{-30,68}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Line(
            points={{16,40},{32,18},{36,-2},{34,-20},{20,-42}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{38,8},{92,-6}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Backward,
            textString="alpha")}));

end SharpEdgedOrifice;

model SuddenExpansion
    "Pressure drop in pipe due to suddenly expanding area (for both flow directions)"
  extends BaseClasses.QuadraticTurbulent.BaseModel(final data=
          BaseClasses.QuadraticTurbulent.LossFactorData.suddenExpansion(
          diameter_a, diameter_b));
  parameter SI.Diameter diameter_a "Inner diameter of pipe at port_a";
  parameter SI.Diameter diameter_b "Inner diameter of pipe at port_b";

  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}},
          grid={1,1}), graphics={
          Line(points={{0,40},{-100,40},{-100,-40},{0,-40},{0,-100},{100,-100},
                {100,100},{0,100},{0,40}}, color={0,0,0}),
          Rectangle(
            extent={{-100,40},{0,-40}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{0,100},{100,-100}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(points={{0,40},{-100,40},{-100,-40},{0,-40},{0,-100},{100,-100},
                {100,100},{0,100},{0,40}}, color={0,0,0}),
          Line(
            points={{-60,-40},{-60,40}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{-50,16},{-26,-10}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="diameter_a"),
          Line(
            points={{34,-100},{34,100}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{54,16},{78,-10}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="diameter_b")}),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}), graphics={
          Text(
            extent={{-150,90},{150,130}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Rectangle(
            extent={{-100,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,20},{0,-20}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Rectangle(
            extent={{0,50},{100,-50}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255})}),
      Documentation(info="<html>
 
</html>"));
end SuddenExpansion;

  model StaticHead
    "Models the static head between two ports at different heights"
    extends Modelica_Fluid.Fittings.BaseClasses.PartialTwoPortTransport;
    outer Modelica_Fluid.System system "System properties";

    parameter SI.Length height_ab "Height(port_b) - Height(port_a)";

    parameter SI.AbsolutePressure dp_nominal=1
      "Nominal pressure drop (0 to disable)" 
      annotation(Dialog(tab="Advanced",group="Regularization"));
    parameter SI.MassFlowRate m_flow_nominal=1 "Mass flow rate for dp_nominal" 
      annotation(Dialog(tab="Advanced",group="Regularization"));

    parameter Boolean smoothFlowReversal=false
      "=true for numerical regularization around zero flow" 
      annotation(Dialog(tab="Advanced",group="Regularization",enable=allowFlowReversal and not use_d_nominal));
    parameter SI.MassFlowRate m_flow_small = 0.01
      "Within regularization if |m_flow| < m_flow_small" 
      annotation(Dialog(tab="Advanced",group="Regularization",enable=smoothFlowReversal));

    parameter Boolean use_d_nominal = false
      "= true to use d_nominal for static head, otherwise computed from medium"
      annotation(Dialog(tab="Advanced"), Evaluate=true);
    parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
      annotation(Dialog(tab="Advanced", enable=use_d_nominal));

  protected
    SI.Density d_a = if use_d_nominal then d_nominal else port_a_d_inflow;
    SI.Density d_b = if use_d_nominal then d_nominal else port_b_d_inflow;
    Real g_times_height_ab(final unit="m2/s2") = system.g*height_ab
      "Gravitiy times height_ab = dp_grav/d";

  equation
    if not allowFlowReversal or use_d_nominal then
      dp = g_times_height_ab*d_a + dp_nominal/m_flow_nominal*m_flow;
    else
      if not smoothFlowReversal then
        // exact switching
        dp = g_times_height_ab*(if m_flow > 0 then d_a else d_b) + dp_nominal/m_flow_nominal*m_flow;
      else
        // regularization around zero flow
        dp = g_times_height_ab*Utilities.regStep(m_flow, d_a, d_b, m_flow_small) + dp_nominal/m_flow_nominal*m_flow;
      end if;
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Rectangle(
            extent={{-100,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,48},{100,-48}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Text(
            extent={{-150,80},{150,120}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name")}),           Documentation(info="<html>
<p>
This model describes the static head due to the relative height between the two connectors. No mass, energy and momentum storage, and no pressure drop due to friction are considered.
</p>
</html>", revisions="<html>
<ul>
<li><i>8 Dec 2008</i>
    by Ruediger Franke:<br>
       Introduce small nominal pressure loss for regularization</li>
<li><i>31 Oct 2007</i>
    by <a href=\"mailto:jonas@modelon.se\">Jonas Eborn</a>:<br>
       Changed to flow-direction dependent density</li>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added to Modelica_Fluid</li>
</ul>
</html>"));
  end StaticHead;

  model WallFrictionAndGravity
    "Pressure drop in pipe due to wall friction and gravity (for both flow directions)"
    extends Modelica_Fluid.Fittings.BaseClasses.PartialTwoPortTransport;

    replaceable package WallFriction = 
      Modelica_Fluid.Fittings.BaseClasses.WallFriction.QuadraticTurbulent 
      constrainedby
      Modelica_Fluid.Fittings.BaseClasses.WallFriction.PartialWallFriction
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
    parameter SI.DynamicViscosity eta_nominal = Medium.dynamicViscosity(
                                                   Medium.setState_pTX(
                                                       Medium.p_default, Medium.T_default, Medium.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                              annotation(Dialog(enable=use_nominal));
    parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(enable=use_nominal));

    parameter Boolean show_Re = false
      "= true, if Reynolds number is included for plotting" 
       annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter Boolean from_dp=true
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
      annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter SI.AbsolutePressure dp_small = 1
      "Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)"
      annotation(Dialog(tab="Advanced", enable=from_dp and WallFriction.use_dp_small));
    parameter SI.MassFlowRate m_flow_small = reg_m_flow_small
      "Within regularizatio if |m_flow| < m_flow_small (may be wider for large discontinuities in static head)"
      annotation(Dialog(tab="Advanced", enable=not from_dp and WallFriction.use_m_flow_small));
    SI.ReynoldsNumber Re = Utilities.ReynoldsNumber_m_flow(m_flow, noEvent(if m_flow>0 then eta_a else eta_b), diameter) if show_Re
      "Reynolds number of pipe";

    outer Modelica_Fluid.System system "System properties";

  protected
    SI.DynamicViscosity eta_a = if not WallFriction.use_eta then 1.e-10 else 
                                (if use_nominal then eta_nominal else Medium.dynamicViscosity(port_a_state_inflow));
    SI.DynamicViscosity eta_b = if not WallFriction.use_eta then 1.e-10 else 
                                (if use_nominal then eta_nominal else Medium.dynamicViscosity(port_b_state_inflow));
    SI.Density d_a = if use_nominal then d_nominal else port_a_d_inflow;
    SI.Density d_b = if use_nominal then d_nominal else port_b_d_inflow;

    Real g_times_height_ab(final unit="m2/s2") = system.g*height_ab
      "Gravitiy times height_ab = dp_grav/d";

    // Currently not in use (means to widen the regularization domain in case of large difference in static head)
    final parameter Boolean use_x_small_staticHead = false
      "Use dp_/m_flow_small_staticHead only if static head actually exists" annotation(Evaluate=true);
                                                           /*abs(height_ab)>0*/
    SI.AbsolutePressure dp_small_staticHead = noEvent(max(dp_small, 0.015*abs(g_times_height_ab*(d_a-d_b))))
      "Heuristic for large discontinuities in static head";
    SI.MassFlowRate m_flow_small_staticHead = noEvent(max(m_flow_small, (-5.55e-7*(d_a+d_b)/2+5.5e-4)*abs(g_times_height_ab*(d_a-d_b))))
      "Heuristic for large discontinuities in static head";

  equation
    if from_dp and not WallFriction.dp_is_zero then
      m_flow = WallFriction.massFlowRate_dp_staticHead(dp, d_a, d_b, eta_a, eta_b, length, diameter,
        g_times_height_ab, roughness, if use_x_small_staticHead then dp_small_staticHead else dp_small);
    else
      dp = WallFriction.pressureLoss_m_flow_staticHead(m_flow, d_a, d_b, eta_a, eta_b, length, diameter,
        g_times_height_ab, roughness, if use_x_small_staticHead then m_flow_small_staticHead else m_flow_small);
    end if;

      annotation (defaultComponentName="pipeFriction",Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,44},{100,-45}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Text(
            extent={{-150,80},{150,120}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name")}),           Documentation(info="<html>
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
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,64},{100,-64}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Backward),
          Rectangle(
            extent={{-100,50},{100,-49}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-60,-49},{-60,50}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{-50,16},{6,-10}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="diameter"),
          Line(
            points={{-100,74},{100,74}},
            color={0,0,255},
            arrow={Arrow.Filled,Arrow.Filled}),
          Text(
            extent={{-34,92},{34,74}},
            lineColor={0,0,255},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="length")}));
  end WallFrictionAndGravity;

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

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

    partial model PartialTJunction
      "Base class for a splitting/joining component with three ports"
      import Modelica_Fluid.Types;
      import Modelica_Fluid.Types.PortFlowDirection;

      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
        "Medium in the component" 
        annotation (choicesAllMatching=true);

      Modelica_Fluid.Interfaces.FluidPort_a port_1(redeclare package Medium = 
            Medium, m_flow(min=if (portFlowDirection_1 == PortFlowDirection.Entering) then 
                    0.0 else -Modelica.Constants.inf, max=if (portFlowDirection_1
               == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
              rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
            Medium, m_flow(min=if (portFlowDirection_2 == PortFlowDirection.Entering) then 
                    0.0 else -Modelica.Constants.inf, max=if (portFlowDirection_2
               == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
        annotation (Placement(transformation(extent={{90,-10},{110,10}}, rotation=
               0)));
      Modelica_Fluid.Interfaces.FluidPort_a port_3(
        redeclare package Medium=Medium,
        m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
        max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
        annotation (Placement(transformation(extent={{-10,90},{10,110}}, rotation=
               0)));

      annotation(Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{-100,41},{100,-47}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={192,192,192}),
            Rectangle(
              extent={{-100,37},{100,-43}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,127,255}),
            Rectangle(
              extent={{-34,100},{34,37}},
              lineColor={0,0,0},
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={192,192,192}),
            Rectangle(
              extent={{-30,100},{30,35}},
              lineColor={0,0,0},
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={0,127,255}),
            Text(
              extent={{-150,-60},{150,-100}},
              lineColor={0,0,255},
              textString="%name")}),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics));

    protected
      parameter PortFlowDirection portFlowDirection_1=PortFlowDirection.Bidirectional
        "Flow direction for port_1" 
       annotation(Dialog(tab="Advanced"));
      parameter PortFlowDirection portFlowDirection_2=PortFlowDirection.Bidirectional
        "Flow direction for port_2" 
       annotation(Dialog(tab="Advanced"));
      parameter PortFlowDirection portFlowDirection_3=PortFlowDirection.Bidirectional
        "Flow direction for port_3" 
       annotation(Dialog(tab="Advanced"));

    end PartialTJunction;

    package SimpleGenericOrifice
      "Simple pressure loss component defined by two constants (diameter, zeta) for the quadratic turbulent regime"

      function massFlowRate_dp
        "Return mass flow rate from pressure drop (m_flow = f(dp))"
        extends Modelica.Icons.Function;
        input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input SI.Diameter D "Diameter at port_a or port_b";
        input Real zeta
          "Constant pressure loss factor with respect to D (i.e., either port_a or port_b)";
        input SI.AbsolutePressure dp_small = 1
          "Turbulent flow if |dp| >= dp_small";
        output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";

        annotation (Documentation(info="<html>
<p>
Compute mass flow rate from constant loss factor and pressure drop (m_flow = f(dp)).
For small pressure drops (dp &lt; dp_small), the characteristic is approximated by 
a polynomial in order to have a finite derivative at zero mass flow rate.
</p>
</html>"));
      algorithm
        /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2 
  */
        m_flow := Utilities.regRoot2(
            dp,
            dp_small,
            d_a/lossConstant_D_zeta(D, zeta),
            d_b/lossConstant_D_zeta(D, zeta));
      end massFlowRate_dp;

      function pressureLoss_m_flow
        "Return pressure drop from mass flow rate (dp = f(m_flow))"
              extends Modelica.Icons.Function;

        input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input SI.Diameter D "Diameter at port_a or port_b";
        input Real zeta
          "Constant pressure loss factor with respect to D (i.e., either port_a or port_b)";
        input SI.MassFlowRate m_flow_small = 0.01
          "Turbulent flow if |m_flow| >= m_flow_small";
        output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";

        annotation (Documentation(info="<html>
<p>
Compute pressure drop from mass flow rate (dp = f(m_flow)).
For small mass flow rates(|m_flow| &lt; m_flow_small), the characteristic is approximated by 
a polynomial in order to have a finite derivative at zero mass flow rate.
</p>
</html>"));
      algorithm
        /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
        dp := Utilities.regSquare2(
            m_flow,
            m_flow_small,
            lossConstant_D_zeta(D, zeta)/d_a,
            lossConstant_D_zeta(D, zeta)/d_b);
      end pressureLoss_m_flow;
      annotation (Documentation(info="<html>
<p>
This pressure drop component defines a
simple, generic orifice, where the loss factor &zeta;=zeta is provided
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
</html>"));
    end SimpleGenericOrifice;

    package QuadraticTurbulent
      "Pressure loss components that are mainly defined by a quadratic turbulent regime with constant loss factor data"
     record LossFactorData
        "Data structure defining constant loss factor data for dp = zeta*rho*v*|v|/2 and functions providing the data for some loss types"

            extends Modelica.Icons.Record;

      SI.Diameter diameter_a "Diameter at port_a" annotation(Dialog);
      SI.Diameter diameter_b "Diameter at port_b" annotation(Dialog);
      Real zeta1 "Loss factor for flow port_a -> port_b" annotation(Dialog);
      Real zeta2 "Loss factor for flow port_b -> port_a" annotation(Dialog);
      SI.ReynoldsNumber Re_turbulent
          "Loss factors suited for Re >= Re_turbulent"                            annotation(Dialog);
      SI.Diameter D_Re "Diameter used to compute Re" annotation(Dialog);
      Boolean zeta1_at_a = true
          "dp = zeta1*(if zeta1_at_a then d_a*v_a^2/2 else d_b*v_b^2/2)" 
                                                                        annotation(Dialog);
      Boolean zeta2_at_a = false
          "dp = -zeta2*(if zeta2_at_a then d_a*v_a^2/2 else d_b*v_b^2/2)" 
                                                                         annotation(Dialog);
      Boolean zetaLaminarKnown = false
          "= true, if zeta = c0/Re in laminar region"                              annotation(Dialog);
      Real c0 = 1
          "zeta = c0/Re; dp = zeta*d_Re*v_Re^2/2, Re=v_Re*D_Re*d_Re/eta_Re)"         annotation(Dialog(enable=zetaLaminarKnown));

      annotation (preferedView="info", Documentation(info="<html>
<p>
This record defines the pressure loss factors of a pipe
segment (orifice, bending etc.) with a minimum amount of data.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
      = 8*&zeta;/(&pi;^2*D^4*&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));

       encapsulated function wallFriction
          "Return pressure loss data due to friction in a straight pipe with walls of nonuniform roughness (not useful for smooth pipes, since zeta is no function of Re)"
          import
            Modelica_Fluid.Fittings.BaseClasses.QuadraticTurbulent.LossFactorData;
          import lg = Modelica.Math.log10;
          import SI = Modelica.SIunits;

         input SI.Length length "Length of pipe" annotation(Dialog);
         input SI.Diameter diameter "Inner diameter of pipe" annotation(Dialog);
         input SI.Length roughness(min=1e-10)
            "Absolute roughness of pipe (> 0 required, details see info layer)"
                                                                               annotation(Dialog);
         output LossFactorData data
            "Pressure loss factors for both flow directions";
         annotation (Icon(coordinateSystem(
                preserveAspectRatio=false,
                extent={{-100,-100},{100,100}},
                grid={1,1}), graphics={Rectangle(
                  extent={{-100,50},{100,-50}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}),
                                   Diagram(coordinateSystem(
                preserveAspectRatio=false,
                extent={{-100,-100},{100,100}},
                grid={1,1}), graphics={
                Rectangle(
                  extent={{-100,64},{100,-64}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Backward),
                Rectangle(
                  extent={{-100,50},{100,-49}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Line(
                  points={{-60,-49},{-60,50}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{-50,16},{6,-10}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="diameter"),
                Line(
                  points={{-100,74},{100,74}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{-34,92},{34,74}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="length")}),
           Documentation(info="<html>
<p>
Friction in straight pipe with walls of nonuniform roughness 
(commercial pipes) in the region that does not depend on the Reynolds-number
</p>
<p>
The loss factors are given for mass flow rates from 
port_a to port_b as:
</p>
<pre>
  turbulent flow (Idelchik 1994, diagram 2-5, p. 117)
     zeta = (L/D)/(2*lg(3.7 / &Delta;))^2, for Re >= 560/&Delta;
&nbsp;
     for Re &ge; 560/&Delta; the loss factor does not depend on the
     Reynolds number. For Re &ge; 4000, the flow is turbulent,
     but depends both on &Delta; and slightly on Re.
&nbsp;
  laminar flow (Idelchick 1994, diagram 2-1, p. 110):
     zeta = 64*(L/D)/Re
</pre>
<p>
where
</p>
<ul>
<li> D is the inner pipe diameter</li>
<li> L is the lenght of the pipe</li>
<li> &Delta; = &delta;/D is the relative roughness where &delta; is
     the absolute \"roughness\", i.e., the averaged height of asperities in the pipe.
     (&delta; may change over time due to growth of surface asperities during
      service, see [Idelchick 1994, p. 85, Tables 2-1, 2-2]).</li>
</ul>
 
<p>
Since the LossFactorData record can only describe loss factors that depend
on geometry (but, e.g., not on the Reynolds number), only the region
with Re &ge; 560/&Delta; is described by this data. Still, the turbulent
region with the above zeta is defined to start at Re=4000, since otherwise
the approximation for Re &lt; 560/&Delta; is too bad.
</p>
 
<p>
The absolute roughness <font face=\"Symbol\">d</font> has usually to
be estimated. In <i>[Idelchik 1994, pp. 105-109,
Table 2-5; Miller 1990, p. 190, Table 8-1]</i> many examples are given.
As a short summary:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Smooth pipes</b></td>
      <td>Drawn brass, coper, aluminium, glass, etc.</td>
      <td><font face=\"Symbol\">d</font> = 0.0025 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Steel pipes</b></td>
      <td>New smooth pipes</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Mortar lined, average finish</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Heavy rust</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Concrete pipes</b></td>
      <td>Steel forms, first class workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Steel forms, average workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Block linings</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
</table>
</html>"));
        protected
         Real Delta = roughness/diameter "relative roughness";
       algorithm
         data.diameter_a          := diameter;
         data.diameter_b          := diameter;
         data.zeta1        := (length/diameter)/(2*lg(3.7 /Delta))^2;
         data.zeta2        := data.zeta1;
         data.Re_turbulent := 4000
            ">= 560/Delta flow does not depend on Re, but interpolation is bad";
         data.D_Re         := diameter;
         data.zeta1_at_a   := true;
         data.zeta2_at_a   := false;
         data.zetaLaminarKnown := true;
         data.c0               := 64*(length/diameter);
       end wallFriction;

       encapsulated function suddenExpansion
          "Return pressure loss data for sudden expansion or contraction in a pipe (for both flow directions)"
          import
            Modelica_Fluid.Fittings.BaseClasses.QuadraticTurbulent.LossFactorData;
          import SI = Modelica.SIunits;
         input SI.Diameter diameter_a "Inner diameter of pipe at port_a" annotation(Dialog);
         input SI.Diameter diameter_b "Inner diameter of pipe at port_b" annotation(Dialog);
         output LossFactorData data
            "Pressure loss factors for both flow directions";
         annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics={
                Rectangle(
                  extent={{-100,40},{0,-40}},
                  lineColor={255,255,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{0,100},{100,-100}},
                  lineColor={255,255,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Line(points={{0,40},{-100,40},{-100,-40},{0,-40},{0,-100},{
                      100,-100},{100,100},{0,100},{0,40}}, color={0,0,0})}),
                                   Diagram(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={
                Line(points={{0,40},{-100,40},{-100,-40},{0,-40},{0,-100},{
                      100,-100},{100,100},{0,100},{0,40}}, color={0,0,0}),
                Rectangle(
                  extent={{-100,40},{0,-40}},
                  lineColor={255,255,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{0,100},{100,-100}},
                  lineColor={255,255,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Line(points={{0,40},{-100,40},{-100,-40},{0,-40},{0,-100},{
                      100,-100},{100,100},{0,100},{0,40}}, color={0,0,0}),
                Line(
                  points={{-60,-40},{-60,40}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{-50,16},{-26,-10}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="diameter_a"),
                Line(
                  points={{34,-100},{34,100}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{54,16},{78,-10}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="diameter_b")}),
           Documentation(info="<html>
<p>
The loss factors are given for mass flow rates from 
port_a to port_b as:
</p>
<pre>
   A_a &lt; A_b (Idelchik 1994, diagram 4-1, p. 208):
      zeta = dp/(d_a*v_a^2/2)
           = (1 - A_a/A_b)^2 for Re_a &ge; 3.3e3 (turbulent flow)
      zeta = 30/Re           for Re_a &lt; 10    (laminar flow)
&nbsp;
   A_a &gt; A_b (Idelchik 1994, diagram 4-9, p. 216 and diagram 4-10, p. 217)
      zeta = dp/(d_b*v_b^2/2)
           = 0.5*(1 - A_b/A_a)^0.75 for Re_b &ge; 1e4 (turbulent flow)
      zeta = 30/Re                  for Re_a &lt; 10  (laminar flow)
</pre>
</html>"));
        protected
         Real A_rel;
       algorithm
         data.diameter_a          := diameter_a;
         data.diameter_b          := diameter_b;
         data.Re_turbulent := 100;
         data.zetaLaminarKnown := true;
         data.c0 := 30;

         if diameter_a <= diameter_b then
            A_rel :=(diameter_a/diameter_b)^2;
            data.zeta1 :=(1 - A_rel)^2;
            data.zeta2 :=0.5*(1 - A_rel)^0.75;
            data.zeta1_at_a :=true;
            data.zeta2_at_a :=true;
            data.D_Re := diameter_a;
         else
            A_rel :=(diameter_b/diameter_a)^2;
            data.zeta1 :=0.5*(1 - A_rel)^0.75;
            data.zeta2 :=(1 - A_rel)^2;
            data.zeta1_at_a :=false;
            data.zeta2_at_a :=false;
            data.D_Re := diameter_b;
         end if;
       end suddenExpansion;

       encapsulated function sharpEdgedOrifice
          "Return pressure loss data for sharp edged orifice (for both flow directions)"
          import NonSI = Modelica.SIunits.Conversions.NonSIunits;
          import
            Modelica_Fluid.Fittings.BaseClasses.QuadraticTurbulent.LossFactorData;
          import SI = Modelica.SIunits;
          input SI.Diameter diameter
            "Inner diameter of pipe (= same at port_a and port_b)" 
                                                                  annotation(Dialog);
          input SI.Diameter leastDiameter "Smallest diameter of orifice" 
                                                                annotation(Dialog);
          input SI.Diameter length "Length of orifice" 
                                                 annotation(Dialog);
          input NonSI.Angle_deg alpha "Angle of orifice" 
                                                        annotation(Dialog);
          output LossFactorData data
            "Pressure loss factors for both flow directions";
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics={
                Rectangle(
                  extent={{-100,60},{100,-60}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{-30,60},{-30,12},{30,50},{30,60},{-30,60}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Backward),
                Polygon(
                  points={{-30,-10},{-30,-60},{30,-60},{30,-50},{-30,-10}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Backward)}),
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics={
                Rectangle(
                  extent={{-100,60},{100,-60}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{-30,60},{-30,12},{30,50},{30,60},{-30,60}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Backward),
                Polygon(
                  points={{-30,-10},{-30,-60},{30,-60},{30,-50},{-30,-10}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Backward),
                Line(
                  points={{-82,-60},{-82,60}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{-78,16},{-44,-8}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="diameter"),
                Line(
                  points={{-30,-10},{-30,12}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{-24,14},{8,-10}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="leastDiameter"),
                Text(
                  extent={{-20,84},{18,70}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid,
                  textString="length"),
                Line(
                  points={{30,68},{-30,68}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Line(
                  points={{16,40},{32,18},{36,-2},{34,-20},{20,-42}},
                  color={0,0,255},
                  arrow={Arrow.Filled,Arrow.Filled}),
                Text(
                  extent={{38,8},{92,-6}},
                  lineColor={0,0,255},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Backward,
                  textString="alpha")}),
            Documentation(info="<html>
<p>
Loss factor for mass flow rate from port_a to port_b
(Idelchik 1994, diagram 4-14, p. 221):
</p>
<pre>
   zeta = [(1-A0/A1) + 0.707*(1-A0/A1)^0.375]^2*(A1/A0)^2 
          for Re(A0) >= 1e5,  independent of alpha
</pre>
<p>
Loss factor for mass flow rate from port_b to port_a
(Idelchik 1994, diagram 4-13, p. 220, with A2=A1):
</p>
<pre>
   zeta = k*(1 - A0/A1)^0.75 + (1 - A0/A1)^2 + 2*sqrt(k*(1-A0/A1)^0.375) + (1- A0/A1)
          k  = 0.13 + 0.34*10^(-(3.4*LD+88.4*LD^2.3)) 
               (there is a typing error in the formula in diagram 4-13, the above
                equation corresponds to table (a) in diagram 4-12)
          LD = L/D0
          for Re(A0) >= 1e4, 40 deg &le; alpha &le; 60 deg 
                             for other values of alpha, k is given as table
                             in diagram 3-7 (this is not yet included in the function)
</pre
</html>"));
        protected
          Real D_rel=leastDiameter/diameter;
          Real LD=length/leastDiameter;
          Real k=0.13 + 0.34*10^(-(3.4*LD + 88.4*LD^2.3));
       algorithm
          data.diameter_a := diameter;
          data.diameter_b := diameter;
          data.zeta1 := ((1 - D_rel) + 0.707*(1 - D_rel)^0.375)^2*(1/D_rel)^2;
          data.zeta2 := k*(1 - D_rel)^0.75 + (1 - D_rel)^2 + 2*sqrt(k*(1 -
            D_rel)^0.375) + (1 - D_rel);
          data.Re_turbulent := 1e4;
          data.D_Re := leastDiameter;
          data.zeta1_at_a := true;
          data.zeta2_at_a := false;
          data.zetaLaminarKnown := false;
          data.c0 := 0;
       end sharpEdgedOrifice;

     end LossFactorData;

      function massFlowRate_dp
        "Return mass flow rate from constant loss factor data and pressure drop (m_flow = f(dp))"
              //import Modelica_Fluid.PressureLosses.BaseClasses.lossConstant_D_zeta;
        extends Modelica.Icons.Function;

        input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input LossFactorData data
          "Constant loss factors for both flow directions" annotation (
            choices(
            choice=Modelica_Fluid.Fittings.Utilities.QuadraticTurbulent.LossFactorData.wallFriction(),
            choice=Modelica_Fluid.Fittings.Utilities.QuadraticTurbulent.LossFactorData.suddenExpansion(),
            choice=Modelica_Fluid.Fittings.Utilities.QuadraticTurbulent.LossFactorData.sharpEdgedOrifice()));
        input SI.AbsolutePressure dp_small = 1
          "Turbulent flow if |dp| >= dp_small";
        output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";

        annotation (smoothOrder=1, Documentation(info="<html>
<p>
Compute mass flow rate from constant loss factor and pressure drop (m_flow = f(dp)).
For small pressure drops (dp &lt; dp_small), the characteristic is approximated by 
a polynomial in order to have a finite derivative at zero mass flow rate.
</p>
</html>"));
      protected
        Real k1 = lossConstant_D_zeta(if data.zeta1_at_a then data.diameter_a else data.diameter_b,data.zeta1);
        Real k2 = lossConstant_D_zeta(if data.zeta2_at_a then data.diameter_a else data.diameter_b,data.zeta2);
      algorithm
        /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
        m_flow :=Utilities.regRoot2(dp, dp_small, d_a/k1, d_b/k2);
      end massFlowRate_dp;
      annotation (Documentation(info="<html>
<p>
This library provides pressure loss factors of a pipe
segment (orifice, bending etc.) with a minimum amount of data.
If available, data can be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
Basically, the pressure drop is defined by the following
equation:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
      = 8*&zeta;/(&pi;^2*D^4*&rho;) * m_flow*|m_flow|
</pre>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<li>
<li> D is the diameter of the pipe segment. If this is not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
</ul>
 
</html>"));

      function massFlowRate_dp_and_Re
        "Return mass flow rate from constant loss factor data, pressure drop and Re (m_flow = f(dp))"
              extends Modelica.Icons.Function;

        input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input SI.DynamicViscosity eta_a "Dynamic viscosity at port_a";
        input SI.DynamicViscosity eta_b "Dynamic viscosity at port_b";
        input LossFactorData data
          "Constant loss factors for both flow directions" annotation (
            choices(
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.wallFriction(),
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.suddenExpansion(),
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.sharpEdgedOrifice()));
        output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";

        annotation (smoothOrder=1, Documentation(info="<html>
<p>
Compute mass flow rate from constant loss factor and pressure drop (m_flow = f(dp)).
If the Reynolds-number Re &ge; data.Re_turbulent, the flow
is treated as a turbulent flow with constant loss factor zeta.
If the Reynolds-number Re &lt; data.Re_turbulent, the flow
is laminar and/or in a transition region between laminar and
turbulent. This region is approximated by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. 
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"data.c0/Re\". 
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
      protected
        constant Real pi=Modelica.Constants.pi;
        Real k0=2*data.c0/(pi*data.D_Re^3);
        Real k1 = lossConstant_D_zeta(if data.zeta1_at_a then data.diameter_a else data.diameter_b,data.zeta1);
        Real k2 = lossConstant_D_zeta(if data.zeta2_at_a then data.diameter_a else data.diameter_b,data.zeta2);
        Real yd0
          "Derivative of m_flow=m_flow(dp) at zero, if data.zetaLaminarKnown";
        SI.AbsolutePressure dp_turbulent
          "The turbulent region is: |dp| >= dp_turbulent";
      algorithm
      /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If data.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
         dp_turbulent :=(k1 + k2)/(d_a + d_b)*
                        ((eta_a + eta_b)*data.D_Re*pi/8)^2*data.Re_turbulent^2;
         yd0 :=if data.zetaLaminarKnown then 
                  (d_a + d_b)/(k0*(eta_a + eta_b)) else 0;
         m_flow := Utilities.regRoot2(dp, dp_turbulent, d_a/k1, d_b/k2,
                                                     data.zetaLaminarKnown, yd0);
      end massFlowRate_dp_and_Re;

      function pressureLoss_m_flow
        "Return pressure drop from constant loss factor and mass flow rate (dp = f(m_flow))"
              extends Modelica.Icons.Function;

        input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input LossFactorData data
          "Constant loss factors for both flow directions" annotation (
            choices(
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.wallFriction(),
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.suddenExpansion(),
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.sharpEdgedOrifice()));
        input SI.MassFlowRate m_flow_small = 0.01
          "Turbulent flow if |m_flow| >= m_flow_small";
        output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";

        annotation (smoothOrder=1, Documentation(info="<html>
<p>
Compute pressure drop from constant loss factor and mass flow rate (dp = f(m_flow)).
For small mass flow rates(|m_flow| &lt; m_flow_small), the characteristic is approximated by 
a polynomial in order to have a finite derivative at zero mass flow rate.
</p>
</html>"));
      protected
        Real k1 = lossConstant_D_zeta(if data.zeta1_at_a then data.diameter_a else data.diameter_b,data.zeta1);
        Real k2 = lossConstant_D_zeta(if data.zeta2_at_a then data.diameter_a else data.diameter_b,data.zeta2);
      algorithm
        /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
        dp :=Utilities.regSquare2(m_flow, m_flow_small, k1/d_a, k2/d_b);
      end pressureLoss_m_flow;

      function pressureLoss_m_flow_and_Re
        "Return pressure drop from constant loss factor, mass flow rate and Re (dp = f(m_flow))"
              extends Modelica.Icons.Function;

        input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input SI.DynamicViscosity eta_a "Dynamic viscosity at port_a";
        input SI.DynamicViscosity eta_b "Dynamic viscosity at port_b";
        input LossFactorData data
          "Constant loss factors for both flow directions" annotation (
            choices(
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.wallFriction(),
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.suddenExpansion(),
            choice=BaseClasses.PressureLosses.QuadraticTurbulent.LossFactorData.sharpEdgedOrifice()));
        output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";

        annotation (smoothOrder=1, Documentation(info="<html>
<p>
Compute pressure drop from constant loss factor and mass flow rate (dp = f(m_flow)).
If the Reynolds-number Re &ge; data.Re_turbulent, the flow
is treated as a turbulent flow with constant loss factor zeta.
If the Reynolds-number Re &lt; data.Re_turbulent, the flow
is laminar and/or in a transition region between laminar and
turbulent. This region is approximated by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. 
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"data.c0/Re\". 
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
      protected
        constant Real pi=Modelica.Constants.pi;
        Real k0 = 2*data.c0/(pi*data.D_Re^3);
        Real k1 = lossConstant_D_zeta(if data.zeta1_at_a then data.diameter_a else data.diameter_b,data.zeta1);
        Real k2 = lossConstant_D_zeta(if data.zeta2_at_a then data.diameter_a else data.diameter_b,data.zeta2);
        Real yd0
          "Derivative of dp = f(m_flow) at zero, if data.zetaLaminarKnown";
        SI.MassFlowRate m_flow_turbulent
          "The turbulent region is: |m_flow| >= m_flow_turbulent";
      algorithm
      /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If data.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
        m_flow_turbulent :=(pi/8)*data.D_Re*(eta_a + eta_b)*data.Re_turbulent;
        yd0 :=if data.zetaLaminarKnown then k0*(eta_a + eta_b)/(d_a + d_b) else 0;
        dp :=Utilities.regSquare2(m_flow, m_flow_turbulent, k1/d_a, k2/d_b,
                                                 data.zetaLaminarKnown, yd0);
      end pressureLoss_m_flow_and_Re;

      partial model BaseModel
        "Generic pressure drop component with constant turbulent loss factor data and without an icon"

        extends Modelica_Fluid.Fittings.BaseClasses.PartialTwoPortTransport;

        SI.ReynoldsNumber Re = Utilities.ReynoldsNumber_m_flow(
              m_flow, (Medium.dynamicViscosity(port_a_state_inflow) + Medium.dynamicViscosity(port_b_state_inflow))/2,
              data.D_Re) if show_Re "Reynolds number at diameter data.D_Re";
        parameter LossFactorData data "Loss factor data";
        parameter Boolean show_Re = false
          "= true, if Reynolds number is included for plotting" 
           annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter Boolean from_dp = true
          "= true, use m_flow = f(dp) else dp = f(m_flow)" 
          annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter Boolean use_Re = false
          "= true, if turbulent region is defined by Re, otherwise by dp_small or m_flow_small"
          annotation(Evaluate=true, Dialog(tab="Advanced"));
        parameter SI.AbsolutePressure dp_small = 1
          "Turbulent flow if |dp| >= dp_small" 
          annotation(Dialog(tab="Advanced", enable=not use_Re and from_dp));
        parameter SI.MassFlowRate m_flow_small = 0.01
          "Turbulent flow if |m_flow| >= m_flow_small" 
          annotation(Dialog(tab="Advanced", enable=not use_Re and not from_dp));

        annotation (
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                  graphics),
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
               graphics),
          Documentation(info="<html>
<p>
This model computes the pressure loss of a pipe
segment (orifice, bending etc.) with a minimum amount of data
provided via parameter <b>data</b>.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
      equation
        if from_dp then
           m_flow = if use_Re then 
                       massFlowRate_dp_and_Re(
                          dp, port_a_d_inflow, port_b_d_inflow,
                          Medium.dynamicViscosity(port_a_state_inflow),
                          Medium.dynamicViscosity(port_b_state_inflow),
                          data) else 
                       massFlowRate_dp(dp, port_a_d_inflow, port_b_d_inflow, data, dp_small);
        else
           dp = if use_Re then 
                   pressureLoss_m_flow_and_Re(
                       m_flow, port_a_d_inflow, port_b_d_inflow,
                       Medium.dynamicViscosity(port_a_state_inflow),
                       Medium.dynamicViscosity(port_b_state_inflow),
                       data) else 
                   pressureLoss_m_flow(m_flow, port_a_d_inflow, port_b_d_inflow, data, m_flow_small);
        end if;
      end BaseModel;

    model TestWallFriction
        "Pressure drop in pipe due to wall friction (only for test purposes; if needed use instead Utilities.WallFriction)"
            extends BaseModel(final data=
              LossFactorData.wallFriction(
              length,
              diameter,
              roughness));
      parameter SI.Length length "Length of pipe";
      parameter SI.Diameter diameter "Inner diameter of pipe";
      parameter SI.Length roughness(min=1e-10)
          "Absolute roughness of pipe (> 0 required, details see info layer)";
      annotation (
        Diagram(graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics={
              Text(
                extent={{-150,80},{150,120}},
                lineColor={0,0,0},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={0,127,255},
                textString="%name"),
              Rectangle(
                extent={{-100,60},{100,-60}},
                lineColor={0,0,0},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={192,192,192}),
              Rectangle(
                extent={{-100,34},{100,-36}},
                lineColor={0,0,0},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={0,127,255}),
              Text(
                extent={{-134,-66},{130,-92}},
                lineColor={0,0,0},
                textString="quad. turbulent")}),
          Documentation(info="<html>
 
</html>"));
    end TestWallFriction;
    end QuadraticTurbulent;

    package WallFriction
      "Different variants for pressure drops due to pipe wall friction"
      partial package PartialWallFriction
        "Partial wall friction characteristic (base package of all wall friction characteristics)"

        annotation (Documentation(info="<html>
 
</html>"));

      // Constants to be set in subpackages
        constant Boolean use_eta = true
          "= true, if eta_a/eta_b are used in function, otherwise value is not used";
        constant Boolean use_roughness = true
          "= true, if roughness is used in function, otherwise value is not used";
        constant Boolean use_dp_small = true
          "= true, if dp_small is used in function, otherwise value is not used";
        constant Boolean use_m_flow_small = true
          "= true, if m_flow_small is used in function, otherwise value is not used";
        constant Boolean dp_is_zero = false
          "= true, if no wall friction is present, i.e., dp = 0 (function massFlowRate_dp() cannot be used)";

      // pressure loss characteristic functions
        replaceable partial function massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          extends Icons.ObsoleteFunction;

          input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.AbsolutePressure dp_small = 1
            "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";

          output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        annotation (Documentation(info="<html>
 
</html>"));
        end massFlowRate_dp;

        replaceable partial function massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
          extends Modelica.Icons.Function;

          input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input Real g_times_height_ab
            "Gravity times (Height(port_b) - Height(port_a))";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.AbsolutePressure dp_small=1
            "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";

          output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          annotation (Documentation(info="<html>
 
</html>"));
        end massFlowRate_dp_staticHead;

        replaceable partial function pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          extends Icons.ObsoleteFunction;

          input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.MassFlowRate m_flow_small = 0.01
            "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
          output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";

        annotation (Documentation(info="<html>
 
</html>"));
        end pressureLoss_m_flow;

        replaceable partial function pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
                  extends Modelica.Icons.Function;

          input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input Real g_times_height_ab
            "Gravity times (Height(port_b) - Height(port_a))";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.MassFlowRate m_flow_small = 0.01
            "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
          output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";

        annotation (Documentation(info="<html>
 
</html>"));
        end pressureLoss_m_flow_staticHead;
      end PartialWallFriction;

      annotation (Documentation(info="<html>
<p>
This package provides functions to compute
pressure losses due to <b>wall friction</b> in a pipe.
Every correlation is defined by a package that is derived
by inheritance from the package WallFriction.PartialWallFriction.
The details of the underlying pipe wall friction model are described in the
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
 
</html>"));
      package NoFriction "No pipe wall friction, no static head"

        annotation (Documentation(info="<html>
<p>
This component sets the pressure loss due to wall friction 
to zero, i.e., it allows to switch off pipe wall friction.
</p>
</html>"));

        extends PartialWallFriction(
                  final use_eta = false,
                  final use_roughness = false,
                  final use_dp_small = false,
                  final use_m_flow_small = false,
                  final dp_is_zero = true);

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          assert(false, "function massFlowRate_dp (option: from_dp=true)
cannot be used for WallFriction.NoFriction. Use instead
function pressureLoss_m_flow (option: from_dp=false)");
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          dp := 0;
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));

        algorithm
          assert(false, "function massFlowRate_dp (option: from_dp=true)
cannot be used for WallFriction.NoFriction. Use instead
function pressureLoss_m_flow (option: from_dp=false)");
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));
        /* To include only static head:
protected 
  Real dp_grav_a = g_times_height_ab*d_a 
    "Static head if mass flows in design direction (a to b)";
  Real dp_grav_b = g_times_height_ab*d_b 
    "Static head if mass flows against design direction (b to a)";
*/
        algorithm
        //  dp := Utilities.regStep(m_flow, dp_grav_a, dp_grav_a, m_flow_small);
          dp := 0;
          assert(abs(g_times_height_ab) < Modelica.Constants.small,
           "WallFriction.NoFriction does not consider static head and cannot be used with height_ab<>0!");
        end pressureLoss_m_flow_staticHead;
      end NoFriction;

      package Laminar
        "Pipe wall friction in the laminar regime (linear correlation)"

        annotation (Documentation(info="<html>
<p>
This component defines only the laminar region of wall friction:
dp = k*m_flow, where \"k\" depends on density and dynamic viscosity.
The roughness of the wall does not have an influence on the laminar
flow and therefore argument roughness is ignored.
Since this is a linear relationship, the occuring systems of equations
are usually much simpler (e.g. either linear instead of non-linear).
By using nominal values for density and dynamic viscosity, the 
systems of equations can still further be reduced. 
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the \"light blue curve\" called
<b>Hagen-Poiseuille</b>.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));

        extends PartialWallFriction(
                  final use_eta = true,
                  final use_roughness = false,
                  final use_dp_small = false,
                  final use_m_flow_small = false);

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          m_flow :=dp*Modelica.Constants.pi*diameter^4*(d_a + d_b)/(128*length*(eta_a + eta_b));
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          dp := m_flow*128*length*(eta_a + eta_b)/(Modelica.Constants.pi*diameter^4*(d_a + d_b));
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));
        protected
          Real k0inv = Modelica.Constants.pi*diameter^4/(128*length)
            "Constant factor";

          Real dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          Real dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real dm_flow_ddp_fric_a = k0inv*d_a/eta_a
            "Slope of mass flow rate over dp if flow in design direction (a to b)";
          Real dm_flow_ddp_fric_b = k0inv*d_b/eta_b
            "Slope of mass flow rate over dp if flow against design direction (b to a)";

          Real dp_a=max(dp_grav_a,dp_grav_b)+dp_small
            "Upper end of regularization domain of the m_flow(dp) relation";
          Real dp_b=min(dp_grav_a,dp_grav_b)-dp_small
            "Lower end of regularization domain of the m_flow(dp) relation";

          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;
        algorithm
        /*
  dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
     = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
     = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
     = 2*c0/(pi*D_Re^3) * eta/d * m_flow
     = k0 * eta/d * m_flow
  k0 = 2*c0/(pi*D_Re^3)
*/

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := dm_flow_ddp_fric_a*(dp-dp_grav_a);
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := dm_flow_ddp_fric_b*(dp-dp_grav_b);
          else
            m_flow_a := dm_flow_ddp_fric_a*(dp_a - dp_grav_a);
            m_flow_b := dm_flow_ddp_fric_b*(dp_b - dp_grav_b);

            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        /*
  m_flow := if dp<dp_b then dm_flow_ddp_b*(dp-dp_grav_b) else 
              (if dp>dp_a then dm_flow_ddp_a*(dp-dp_grav_a) else 
                Modelica_Fluid.Utilities.regFun3(dp, dp_b, dp_a, dm_flow_ddp_b*(dp_b - dp_grav_b), dm_flow_ddp_a*(dp_a - dp_grav_a), dm_flow_ddp_b, dm_flow_ddp_a));
*/
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));
        protected
          Real k0 = 128*length/(Modelica.Constants.pi*diameter^4)
            "Constant factor";

          Real dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          Real dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real ddp_dm_flow_a = k0*eta_a/d_a
            "Slope of dp over mass flow rate if flow in design direction (a to b)";
          Real ddp_dm_flow_b = k0*eta_b/d_b
            "Slope of dp over mass flow rate if flow against design direction (b to a)";

          Real m_flow_a=if dp_grav_a >= dp_grav_b then m_flow_small else m_flow_small + (dp_grav_b-dp_grav_a)/ddp_dm_flow_a
            "Upper end of regularization domain of the dp(m_flow) relation";
          Real m_flow_b=if dp_grav_a >= dp_grav_b then -m_flow_small else -m_flow_small - (dp_grav_b - dp_grav_a)/ddp_dm_flow_b
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;
        algorithm
        /*
  dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
     = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
     = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
     = 2*c0/(pi*D_Re^3) * eta/d * m_flow
     = k0 * eta/d * m_flow
  k0 = 2*c0/(pi*D_Re^3)
*/

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := (ddp_dm_flow_a*m_flow + dp_grav_a);
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := (ddp_dm_flow_b*m_flow + dp_grav_b);
          else
            // Regularization parameters
            dp_a := ddp_dm_flow_a*m_flow_a + dp_grav_a;
            dp_b := ddp_dm_flow_b*m_flow_b + dp_grav_b;
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;
        end pressureLoss_m_flow_staticHead;
      end Laminar;

      package QuadraticTurbulent
        "Pipe wall friction in the quadratic turbulent regime (simple characteristic, eta not used)"

        annotation (Documentation(info="<html>
<p>
This component defines only the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is no longer a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the asymptotic behaviour for large
Reynolds numbers, i.e., the values at the right ordinate where
&lambda; is constant.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));

        extends PartialWallFriction(
                  final use_eta = false,
                  final use_roughness = true,
                  final use_dp_small = true,
                  final use_m_flow_small = true);

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta;
          Real k_inv;
        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta  := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k_inv := (pi*diameter*diameter)^2/(8*zeta);
          m_flow := Modelica_Fluid.Utilities.regRoot2(dp, dp_small, d_a*k_inv, d_b*k_inv);
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          import Modelica.Math;

          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta;
          Real k;
        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k    := 8*zeta/(pi*diameter*diameter)^2;
          dp   := Modelica_Fluid.Utilities.regSquare2(m_flow, m_flow_small, k/d_a, k/d_b);
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta = (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          Real k_inv = (pi*diameter*diameter)^2/(8*zeta);

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real k1 = d_a*k_inv "Factor in m_flow =  sqrt(k1*(dp-dp_grav_a))";
          Real k2 = d_b*k_inv "Factor in m_flow = -sqrt(k2*|dp-dp_grav_b|)";

          Real dp_a=max(dp_grav_a,dp_grav_b)+dp_small
            "Upper end of regularization domain of the m_flow(dp) relation";
          Real dp_b=min(dp_grav_a,dp_grav_b)-dp_small
            "Lower end of regularization domain of the m_flow(dp) relation";

          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          SI.MassFlowRate dm_flow_ddp_fric_a
            "Derivative at upper end of regularization domain";
          SI.MassFlowRate dm_flow_ddp_fric_b
            "Derivative at lower end of regularization domain";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;
        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := sqrt(k1*(dp-dp_grav_a));
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := -sqrt(k2*abs(dp-dp_grav_b));
          else
            m_flow_a := sqrt(k1*(dp_a - dp_grav_a));
            m_flow_b := -sqrt(k2*abs(dp_b - dp_grav_b));

            dm_flow_ddp_fric_a := k1/(2*sqrt(k1*(dp_a - dp_grav_a)));
            dm_flow_ddp_fric_b := k2/(2*sqrt(k2*abs(dp_b - dp_grav_b)));
        /*  dm_flow_ddp_fric_a := if abs(dp_a - dp_grav_a)>0 then k1/(2*sqrt(k1*(dp_a - dp_grav_a))) else  Modelica.Constants.inf);
    dm_flow_ddp_fric_b := if abs(dp_b - dp_grav_b)>0 then k2/(2*sqrt(k2*abs(dp_b - dp_grav_b))) else Modelica.Constants.inf; */

            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta = (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          Real k = 8*zeta/(pi*diameter*diameter)^2;

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real k1 = k/d_a "If m_flow >= 0 then dp = k1*m_flow^2 + dp_grav_a";
          Real k2 = k/d_b "If m_flow < 0 then dp = -k2*m_flow^2 + dp_grav_b";

          Real m_flow_a=if dp_grav_a >= dp_grav_b then m_flow_small else m_flow_small + sqrt((dp_grav_b - dp_grav_a)/k1)
            "Upper end of regularization domain of the dp(m_flow) relation";
          Real m_flow_b=if dp_grav_a >= dp_grav_b then -m_flow_small else -m_flow_small - sqrt((dp_grav_b - dp_grav_a)/k2)
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          Real ddp_dm_flow_a
            "Derivative of pressure drop with mass flow rate at m_flow_a";
          Real ddp_dm_flow_b
            "Derivative of pressure drop with mass flow rate at m_flow_b";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;

        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := (k1*m_flow^2 + dp_grav_a);
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := (-k2*m_flow^2 + dp_grav_b);
          else
            // Regularization parameters
            dp_a := k1*m_flow_a^2 + dp_grav_a;
            ddp_dm_flow_a := 2*k1*m_flow_a;
            dp_b := -k2*m_flow_b^2 + dp_grav_b;
            ddp_dm_flow_b := -2*k2*m_flow_b;
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;

        end pressureLoss_m_flow_staticHead;
      end QuadraticTurbulent;

      package LaminarAndQuadraticTurbulent
        "Pipe wall friction in the laminar and quadratic turbulent regime (simple characteristic)"

        annotation (Documentation(info="<html>
<p>
This component defines the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is no longer a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
At Re=4000, a polynomial is constructed that approaches
the constant &lambda; (for large Reynolds-numbers) at Re=4000
smoothly and has a derivative at zero mass flow rate that is
identical to laminar wall friction.
</p>
</html>"));

        extends PartialWallFriction(
                  final use_eta = true,
                  final use_roughness = true,
                  final use_dp_small = true,
                  final use_m_flow_small = true);

        import ln = Modelica.Math.log "Logarithm, base e";
        import Modelica.Math.log10 "Logarithm, base 10";
        import Modelica.Math.exp "Exponential function";
        import Modelica.Constants.pi;

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi=Modelica.Constants.pi;
          constant Real Re_turbulent = 4000 "Start of turbulent regime";
          Real zeta;
          Real k0;
          Real k_inv;
          Real yd0 "Derivative of m_flow=m_flow(dp) at zero";
          SI.AbsolutePressure dp_turbulent;
        algorithm
        /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If data.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta   := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k0     := 128*length/(pi*diameter^4);
          k_inv  := (pi*diameter*diameter)^2/(8*zeta);
          yd0    := (d_a + d_b)/(k0*(eta_a + eta_b));
          dp_turbulent := ((eta_a + eta_b)*diameter*pi/8)^2*Re_turbulent^2/(k_inv*(d_a+d_b)/2);
          m_flow := Modelica_Fluid.Utilities.regRoot2(dp, dp_turbulent, d_a*k_inv, d_b*k_inv,
                                                      use_yd0=true, yd0=yd0);
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          import Modelica.Math;

          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi=Modelica.Constants.pi;
          constant Real Re_turbulent = 4000 "Start of turbulent regime";
          Real zeta;
          Real k0;
          Real k;
          Real yd0 "Derivative of dp = f(m_flow) at zero";
          SI.MassFlowRate m_flow_turbulent
            "The turbulent region is: |m_flow| >= m_flow_turbulent";

        algorithm
        /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
*/
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k0   := 128*length/(pi*diameter^4);
          k    := 8*zeta/(pi*diameter*diameter)^2;
          yd0  := k0*(eta_a + eta_b)/(d_a + d_b);
          m_flow_turbulent :=(pi/8)*diameter*(eta_a + eta_b)*Re_turbulent;
          dp :=Modelica_Fluid.Utilities.regSquare2(m_flow, m_flow_turbulent, k/d_a, k/d_b,
                                                   use_yd0=true, yd0=yd0);
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";

          SI.Pressure dp_a
            "Upper end of regularization domain of the m_flow(dp) relation";
          SI.Pressure dp_b
            "Lower end of regularization domain of the m_flow(dp) relation";

          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          SI.MassFlowRate dm_flow_ddp_fric_a
            "Derivative at upper end of regularization domain";
          SI.MassFlowRate dm_flow_ddp_fric_b
            "Derivative at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;
        algorithm
          assert(roughness > 1.e-10,
            "roughness > 0 required for quadratic turbulent wall friction characteristic");

          dp_a := max(dp_grav_a, dp_grav_b)+dp_small;
          dp_b := min(dp_grav_a, dp_grav_b)-dp_small;

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp - dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          else
            // Regularization parameters
            (m_flow_a, dm_flow_ddp_fric_a) := Internal.m_flow_of_dp_fric(dp_a-dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            (m_flow_b, dm_flow_ddp_fric_b) := Internal.m_flow_of_dp_fric(dp_b-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";

          SI.MassFlowRate m_flow_a
            "Upper end of regularization domain of the dp(m_flow) relation";
          SI.MassFlowRate m_flow_b
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real ddp_dm_flow_a
            "Derivative of pressure drop with mass flow rate at m_flow_a";
          Real ddp_dm_flow_b
            "Derivative of pressure drop with mass flow rate at m_flow_b";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;

        algorithm
          assert(roughness > 1.e-10,
            "roughness > 0 required for quadratic turbulent wall friction characteristic");

          m_flow_a := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_b - dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta)+m_flow_small else 
            m_flow_small;
          m_flow_b := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_a - dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta)-m_flow_small else 
            -m_flow_small;

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_a;
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_b;
          else
            // Regularization parameters
            (dp_a, ddp_dm_flow_a) := Internal.dp_fric_of_m_flow(m_flow_a, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta);
            dp_a := dp_a + dp_grav_a "Adding dp_grav to dp_fric to get dp";
            (dp_b, ddp_dm_flow_b) := Internal.dp_fric_of_m_flow(m_flow_b, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta);
            dp_b := dp_b + dp_grav_b "Adding dp_grav to dp_fric to get dp";
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;
        end pressureLoss_m_flow_staticHead;

        package Internal
          "Functions to calculate mass flow rate from friction pressure drop and vice versa"
          function m_flow_of_dp_fric
            "Calculate mass flow rate as function of pressure drop due to friction"

            input SI.Pressure dp_fric
              "Pressure drop due to friction (dp = port_a.p - port_b.p)";
            input SI.Density d_a "Density at port_a";
            input SI.Density d_b "Density at port_b";
            input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
            input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
            input SI.Length length "Length of pipe";
            input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
            input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            output SI.MassFlowRate m_flow
              "Mass flow rate from port_a to port_b";
            output Real dm_flow_ddp_fric
              "Derivative of mass flow rate with dp_fric";
            annotation (smoothOrder=1);
          protected
            SI.DynamicViscosity eta "Upstream viscosity";
            SI.Density d "Upstream density";

            Real zeta;
            Real k0;
            Real k_inv;
            Real dm_flow_ddp_laminar
              "Derivative of m_flow=m_flow(dp) in laminar regime";
            SI.AbsolutePressure dp_fric_turbulent
              "The turbulent region is: |dp_fric| >= dp_fric_turbulent, simple quadratic correlation";
            SI.AbsolutePressure dp_fric_laminar
              "The laminar region is: |dp_fric| <= dp_fric_laminar";
          algorithm
          /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   dp_fric_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
*/
            // Determine upstream density and upstream viscosity
            if dp_fric >= 0 then
              d := d_a;
              eta := eta_a;
            else
              d := d_b;
              eta := eta_b;
            end if;
            // Quadratic turbulent
            zeta := (length/diameter)/(2*log10(3.7/(Delta)))^2;
            k_inv := (pi*diameter*diameter)^2/(8*zeta);
            dp_fric_turbulent := sign(dp_fric)*(eta*diameter*pi/4)^2*Re2^2/(k_inv*d);

            // Laminar
            k0 := 128*length/(pi*diameter^4);
            dm_flow_ddp_laminar := d/(k0*eta);
            dp_fric_laminar := sign(dp_fric)*pi*k0*eta^2/d*diameter/4*Re1;

            if abs(dp_fric) > abs(dp_fric_turbulent) then
              m_flow := sign(dp_fric)*sqrt(d*k_inv*abs(dp_fric));
              dm_flow_ddp_fric := 0.5*d*k_inv*(d*k_inv*abs(dp_fric))^(-0.5);
            elseif abs(dp_fric) < abs(dp_fric_laminar) then
              m_flow := dm_flow_ddp_laminar*dp_fric;
              dm_flow_ddp_fric := dm_flow_ddp_laminar;
            else
              // Preliminary testing seems to indicate that the log-log transform is not required here
              (m_flow,dm_flow_ddp_fric) := Utilities.cubicHermite_withDerivative(
                dp_fric, dp_fric_laminar, dp_fric_turbulent, dm_flow_ddp_laminar*dp_fric_laminar,
                sign(dp_fric_turbulent)*sqrt(d*k_inv*abs(dp_fric_turbulent)), dm_flow_ddp_laminar,
                if abs(dp_fric_turbulent)>0 then 0.5*d*k_inv*(d*k_inv*abs(dp_fric_turbulent))^(-0.5) else Modelica.Constants.inf);
            end if;
          end m_flow_of_dp_fric;

          function dp_fric_of_m_flow
            "Calculate pressure drop due to friction as function of mass flow rate"

            input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
            input SI.Density d_a "Density at port_a";
            input SI.Density d_b "Density at port_b";
            input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
            input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
            input SI.Length length "Length of pipe";
            input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
            input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            output SI.Pressure dp_fric
              "Pressure drop due to friction (dp_fric = port_a.p - port_b.p - dp_grav)";
            output Real ddp_fric_dm_flow
              "Derivative of pressure drop with mass flow rate";
            annotation (smoothOrder=1);
          protected
            SI.DynamicViscosity eta "Upstream viscosity";
            SI.Density d "Upstream density";
            Real zeta;
            Real k0;
            Real k;
            Real ddp_fric_dm_flow_laminar
              "Derivative of dp_fric = f(m_flow) at zero";
            SI.MassFlowRate m_flow_turbulent
              "The turbulent region is: |m_flow| >= m_flow_turbulent";
            SI.MassFlowRate m_flow_laminar
              "The laminar region is: |m_flow| <= m_flow_laminar";
          algorithm
          /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
*/
            // Determine upstream density and upstream viscosity
            if m_flow >= 0 then
              d := d_a;
              eta := eta_a;
            else
              d := d_b;
              eta := eta_b;
            end if;

            // Turbulent
            zeta := (length/diameter)/(2*log10(3.7/(Delta)))^2;
            k := 8*zeta/(pi*diameter*diameter)^2;
            m_flow_turbulent := sign(m_flow)*(pi/4)*diameter*eta*Re2;

            // Laminar
            k0 := 128*length/(pi*diameter^4);
            ddp_fric_dm_flow_laminar := k0*eta/d;
            m_flow_laminar := sign(m_flow)*(pi/4)*diameter*eta*Re1;

            if abs(m_flow) > abs(m_flow_turbulent) then
              dp_fric := k/d*m_flow*abs(m_flow);
              ddp_fric_dm_flow := 2*k/d*abs(m_flow);
            elseif abs(m_flow) < abs(m_flow_laminar) then
              dp_fric := ddp_fric_dm_flow_laminar*m_flow;
              ddp_fric_dm_flow := ddp_fric_dm_flow_laminar;
            else
              // Preliminary testing seems to indicate that the log-log transform is not required here
              (dp_fric,ddp_fric_dm_flow) := Utilities.cubicHermite_withDerivative(
                m_flow, m_flow_laminar, m_flow_turbulent, ddp_fric_dm_flow_laminar*m_flow_laminar,
                k/d*m_flow_turbulent*abs(m_flow_turbulent), ddp_fric_dm_flow_laminar, 2*k/d*abs(m_flow_turbulent));
            end if;
          end dp_fric_of_m_flow;
        end Internal;
      end LaminarAndQuadraticTurbulent;

      package Detailed
        "Pipe wall friction in the whole regime (detailed characteristic)"

        annotation (Documentation(info="<html>
<p>
This component defines the complete regime of wall friction.
The details are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
The functional relationship of the friction loss factor &lambda; is
displayed in the next figure. Function massFlowRate_dp() defines the \"red curve\"
(\"Swamee and Jain\"), where as function pressureLoss_m_flow() defines the
\"blue curve\" (\"Colebrook-White\"). The two functions are inverses from 
each other and give slightly different results in the transition region
between Re = 1500 .. 4000, in order to get explicit equations without
solving a non-linear equation.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
<p>
Additionally to wall friction, this component properly implements static
head. With respect to the latter, two cases can be distuinguised. In the case
shown next, the change of elevation with the path from a to b has the opposite
sign of the change of density.</p>
 
<img src=\"../Images/Components/PipeFrictionStaticHead_case-a.PNG\">
 
<p>
In the case illustrated second, the change of elevation with the path from a to 
b has the same sign of the change of density.</p>
 
<img src=\"../Images/Components/PipeFrictionStaticHead_case-b.PNG\">
 
 
</html>"));

        extends PartialWallFriction(
                  final use_eta = true,
                  final use_roughness = true,
                  final use_dp_small = true,
                  final use_m_flow_small = true);

        import ln = Modelica.Math.log "Logarithm, base e";
        import Modelica.Math.log10 "Logarithm, base 10";
        import Modelica.Math.exp "Exponential function";
        import Modelica.Constants.pi;

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          import Modelica.Math;
                  annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = (745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97
            "Re leaving laminar curve";
          SI.ReynoldsNumber Re2 = 4000 "Re entering turbulent curve";
          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          SI.ReynoldsNumber Re "Reynolds number";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";

          function interpolateInRegion2
             input Real Re_turbulent;
             input SI.ReynoldsNumber Re1;
             input SI.ReynoldsNumber Re2;
             input Real Delta;
             input Real lambda2;
             output SI.ReynoldsNumber Re;
             annotation(smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
          protected
            Real x1=Math.log10(64*Re1);
            Real y1=Math.log10(Re1);
            Real yd1=1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux1=(0.5/Math.log(10))*5.74*0.9;
            Real aux2=Delta/3.7 + 5.74/Re2^0.9;
            Real aux3=Math.log10(aux2);
            Real L2=0.25*(Re2/aux3)^2;
            Real aux4=2.51/sqrt(L2) + 0.27*Delta;
            Real aux5=-2*sqrt(L2)*Math.log10(aux4);
            Real x2=Math.log10(L2);
            Real y2=Math.log10(aux5);
            Real yd2=0.5 + (2.51/Math.log(10))/(aux5*aux4);

            // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
            Real diff_x=x2 - x1;
            Real m=(y2 - y1)/diff_x;
            Real c2=(3*m - 2*yd1 - yd2)/diff_x;
            Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
            Real lambda2_1=64*Re1;
            Real dx;
          algorithm
             dx := Math.log10(lambda2/lambda2_1);
             Re := Re1*(lambda2/lambda2_1)^(1 + dx*(c2 + dx*c3));
          end interpolateInRegion2;

        algorithm
          // Determine upstream density, upstream viscosity, and lambda2
          d       := if dp >= 0 then d_a else d_b;
          eta     := if dp >= 0 then eta_a else eta_b;
          lambda2 := abs(dp)*2*diameter^3*d/(length*eta*eta);

          // Determine Re under the assumption of laminar flow
          Re := lambda2/64;

          // Modify Re, if turbulent flow
          if Re > Re1 then
             Re :=-2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
             if Re < Re2 then
                Re := interpolateInRegion2(Re, Re1, Re2, Delta, lambda2);
             end if;
          end if;

          // Determine mass flow rate
          m_flow := (pi*diameter/4)*eta*(if dp >= 0 then Re else -Re);
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          import Modelica.Math;
                  annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Re leaving laminar curve";
          SI.ReynoldsNumber Re2 = 4000 "Re entering turbulent curve";
          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          SI.ReynoldsNumber Re "Reynolds number";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";

          function interpolateInRegion2
             input SI.ReynoldsNumber Re;
             input SI.ReynoldsNumber Re1;
             input SI.ReynoldsNumber Re2;
             input Real Delta;
             output Real lambda2;
             annotation(smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
          protected
            Real x1 = Math.log10(Re1);
            Real y1 = Math.log10(64*Re1);
            Real yd1=1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux1=(0.5/Math.log(10))*5.74*0.9;
            Real aux2=Delta/3.7 + 5.74/Re2^0.9;
            Real aux3=Math.log10(aux2);
            Real L2=0.25*(Re2/aux3)^2;
            Real aux4=2.51/sqrt(L2) + 0.27*Delta;
            Real aux5=-2*sqrt(L2)*Math.log10(aux4);
            Real x2 =  Math.log10(Re2);
            Real y2 =  Math.log10(L2);
            Real yd2 = 2 + 4*aux1/(aux2*aux3*(Re2)^0.9);

            // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
            Real diff_x=x2 - x1;
            Real m=(y2 - y1)/diff_x;
            Real c2=(3*m - 2*yd1 - yd2)/diff_x;
            Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
            Real dx;
          algorithm
             dx := Math.log10(Re/Re1);
             lambda2 := 64*Re1*(Re/Re1)^(1 + dx*(c2 + dx*c3));
          end interpolateInRegion2;
        algorithm
          // Determine upstream density and upstream viscosity
          d       :=if m_flow >= 0 then d_a else d_b;
          eta     :=if m_flow >= 0 then eta_a else eta_b;

          // Determine Re, lambda2 and pressure drop
          Re :=(4/pi)*abs(m_flow)/(diameter*eta);
          lambda2 := if Re <= Re1 then 64*Re else 
                    (if Re >= Re2 then 0.25*(Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else 
                     interpolateInRegion2(Re, Re1, Re2, Delta));
          dp :=length*eta*eta/(2*d*diameter*diameter*diameter)*
               (if m_flow >= 0 then lambda2 else -lambda2);
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"

          annotation (smoothOrder=1);

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re "Reynolds number";
          SI.ReynoldsNumber Re1 = (745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";
          SI.Pressure dp_a
            "Upper end of regularization domain of the m_flow(dp) relation";
          SI.Pressure dp_b
            "Lower end of regularization domain of the m_flow(dp) relation";
          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          SI.MassFlowRate dm_flow_ddp_fric_a
            "Derivative at upper end of regularization domain";
          SI.MassFlowRate dm_flow_ddp_fric_b
            "Derivative at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;

        algorithm
          dp_a := max(dp_grav_a, dp_grav_b)+dp_small;
          dp_b := min(dp_grav_a, dp_grav_b)-dp_small;

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          else
            // Regularization parameters
            (m_flow_a, dm_flow_ddp_fric_a) := Internal.m_flow_of_dp_fric(dp_a-dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            (m_flow_b, dm_flow_ddp_fric_b) := Internal.m_flow_of_dp_fric(dp_b-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"

          annotation (smoothOrder=1);

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";

          SI.MassFlowRate m_flow_a
            "Upper end of regularization domain of the dp(m_flow) relation";
          SI.MassFlowRate m_flow_b
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real ddp_dm_flow_a
            "Derivative of pressure drop with mass flow rate at m_flow_a";
          Real ddp_dm_flow_b
            "Derivative of pressure drop with mass flow rate at m_flow_b";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;

        algorithm
          m_flow_a := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_b - dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta)+m_flow_small else 
            m_flow_small;
          m_flow_b := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_a - dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta)-m_flow_small else 
            -m_flow_small;

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_a;
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_b;
          else
            // Regularization parameters
            (dp_a, ddp_dm_flow_a) := Internal.dp_fric_of_m_flow(m_flow_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            dp_a := dp_a + dp_grav_a "Adding dp_grav to dp_fric to get dp";
            (dp_b, ddp_dm_flow_b) := Internal.dp_fric_of_m_flow(m_flow_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            dp_b := dp_b + dp_grav_b "Adding dp_grav to dp_fric to get dp";
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;
        end pressureLoss_m_flow_staticHead;

      package Internal
          "Functions to calculate mass flow rate from friction pressure drop and vice versa"
        function m_flow_of_dp_fric
            "Calculate mass flow rate as function of pressure drop due to friction"

          input SI.Pressure dp_fric
              "Pressure drop due to friction (dp = port_a.p - port_b.p)";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
          input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
          input Real Delta "Relative roughness";
          output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          output Real dm_flow_ddp_fric
              "Derivative of mass flow rate with dp_fric";
          annotation(smoothOrder=1);

          protected
          function interpolateInRegion2_withDerivative
              "Interpolation in log-log space using a cubic Hermite polynomial, where x=log10(lambda2), y=log10(Re)"

            input Real lambda2 "Known independent variable";
            input SI.ReynoldsNumber Re1
                "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
                "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            input SI.Pressure dp_fric
                "Pressure drop due to friction (dp = port_a.p - port_b.p)";
            output SI.ReynoldsNumber Re "Unknown return variable";
            output Real dRe_ddp "Derivative of return value";
            annotation (smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
            protected
            Real x1=log10(64*Re1);
            Real y1=log10(Re1);
            Real y1d=1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux2=Delta/3.7 + 5.74/Re2^0.9;
            Real aux3=log10(aux2);
            Real L2=0.25*(Re2/aux3)^2;
            Real aux4=2.51/sqrt(L2) + 0.27*Delta;
            Real aux5=-2*sqrt(L2)*log10(aux4);
            Real x2=log10(L2);
            Real y2=log10(aux5);
            Real y2d=0.5 + (2.51/ln(10))/(aux5*aux4);

            // Point of interest in transformed space
            Real x=log10(lambda2);
            Real y;
            Real dy_dx "Derivative in transformed space";
          algorithm
            // Interpolation
            (y, dy_dx) := Utilities.cubicHermite_withDerivative(x, x1, x2, y1, y2, y1d, y2d);

            // Return value
            Re := 10^y;

            // Derivative of return value
            dRe_ddp := Re/abs(dp_fric)*dy_dx;
          end interpolateInRegion2_withDerivative;

          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
          SI.ReynoldsNumber Re "Reynolds number";
          Real dRe_ddp "dRe/ddp";
          Real aux1;
          Real aux2;

        algorithm
          // Determine upstream density and upstream viscosity
          if dp_fric >= 0 then
            d := d_a;
            eta := eta_a;
          else
            d := d_b;
            eta := eta_b;
          end if;

          // Positive mass flow rate
          lambda2 := abs(dp_fric)*2*diameter^3*d/(length*eta*eta)
              "Known as lambda2=f(dp)";

          aux1:=(2*diameter^3*d)/(length*eta^2);

          // Determine Re and dRe/ddp under the assumption of laminar flow
          Re := lambda2/64 "Hagen-Poiseuille";
          dRe_ddp := aux1/64 "Hagen-Poiseuille";

          // Modify Re, if turbulent flow
          if Re > Re1 then
            Re :=-2*sqrt(lambda2)*log10(2.51/sqrt(lambda2) + 0.27*Delta)
                "Colebrook-White";
            aux2 := sqrt(aux1*abs(dp_fric));
            dRe_ddp := 1/ln(10)*(-2*ln(2.51/aux2+0.27*Delta)*aux1/(2*aux2)+2*2.51/(2*abs(dp_fric)*(2.51/aux2+0.27*Delta)));
            if Re < Re2 then
              (Re, dRe_ddp) := interpolateInRegion2_withDerivative(lambda2, Re1, Re2, Delta, dp_fric);
            end if;
          end if;

          // Determine mass flow rate
          m_flow := (pi*diameter/4)*eta*(if dp_fric >= 0 then Re else -Re);
          // Determine derivative of mass flow rate with dp_fric
          dm_flow_ddp_fric := (pi*diameter*eta)/4*dRe_ddp;
        end m_flow_of_dp_fric;

        function dp_fric_of_m_flow
            "Calculate pressure drop due to friction as function of mass flow rate"

          input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
          input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
          input Real Delta "Relative roughness";
          output SI.Pressure dp_fric
              "Pressure drop due to friction (dp_fric = port_a.p - port_b.p - dp_grav)";
          output Real ddp_fric_dm_flow
              "Derivative of pressure drop with mass flow rate";
          annotation(smoothOrder=1);

          protected
          function interpolateInRegion2
              "Interpolation in log-log space using a cubic Hermite polynomial, where x=log10(Re), y=log10(lambda2)"

            input SI.ReynoldsNumber Re "Known independent variable";
            input SI.ReynoldsNumber Re1
                "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
                "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
            output Real lambda2 "Unknown return value";
            output Real dlambda2_dm_flow "Derivative of return value";
            annotation(smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
            protected
            Real x1 = log10(Re1);
            Real y1 = log10(64*Re1);
            Real y1d = 1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux2 = Delta/3.7 + 5.74/Re2^0.9;
            Real aux3 = log10(aux2);
            Real L2 = 0.25*(Re2/aux3)^2;
            Real x2 = log10(Re2);
            Real y2 = log10(L2);
            Real y2d = 2+(2*5.74*0.9)/(ln(aux2)*Re2^0.9*aux2);

            // Point of interest in transformed space
            Real x=log10(Re);
            Real y;
            Real dy_dx "Derivative in transformed space";
          algorithm
            // Interpolation
            (y, dy_dx) := Utilities.cubicHermite_withDerivative(x, x1, x2, y1, y2, y1d, y2d);

            // Return value
            lambda2 := 10^y;

            // Derivative of return value
            dlambda2_dm_flow := lambda2/abs(m_flow)*dy_dx;
          end interpolateInRegion2;

          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          SI.ReynoldsNumber Re "Reynolds number";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
          Real dlambda2_dm_flow "dlambda2/dm_flow";
          Real aux1;
          Real aux2;

        algorithm
          // Determine upstream density and upstream viscosity
          if m_flow >= 0 then
            d := d_a;
            eta := eta_a;
          else
            d := d_b;
            eta := eta_b;
          end if;

          // Determine Reynolds number
          Re :=(4/pi)*abs(m_flow)/(diameter*eta);

          aux1 := 4/(pi*diameter*eta);

          // Use correlation for lambda2 depending on actual conditions
          if Re <= Re1 then
            lambda2 := 64*Re "Hagen-Poiseuille";
            dlambda2_dm_flow := 64*aux1 "Hagen-Poiseuille";
          elseif Re >= Re2 then
            lambda2 := 0.25*(Re/log10(Delta/3.7 + 5.74/Re^0.9))^2 "Swamee-Jain";
            aux2 := Delta/3.7+5.74/((aux1*abs(m_flow))^0.9);
            dlambda2_dm_flow := 0.5*aux1*Re*ln(10)^2*(1/(ln(aux2)^2)+(5.74*0.9)/(ln(aux2)^3*Re^0.9*aux2))
                "Swamee-Jain";
          else
            (lambda2, dlambda2_dm_flow) := interpolateInRegion2(Re, Re1, Re2, Delta, m_flow);
          end if;

          // Compute pressure drop from lambda2
          dp_fric :=length*eta*eta/(2*d*diameter*diameter*diameter)*
               (if m_flow >= 0 then lambda2 else -lambda2);

          // Compute derivative from dlambda2/dm_flow
          ddp_fric_dm_flow := (length*eta^2)/(2*diameter^3*d)*dlambda2_dm_flow;
        end dp_fric_of_m_flow;
      end Internal;
      end Detailed;
    end WallFriction;

    function lossConstant_D_zeta "Return the loss constant 8*zeta/(pi^2*D^4)"
          extends Modelica.Icons.Function;

      input SI.Diameter D "Diameter at port_a or port_b";
      input Real zeta
        "Constant pressure loss factor with respect to D (i.e., either port_a or port_b)";
      output Real k "Loss constant (= 8*zeta/(pi^2*D^4))";

      annotation (Documentation(info="<html>
 
</html>"));
    algorithm
      k :=8*zeta/(Modelica.Constants.pi*Modelica.Constants.pi*D*D*D*D);
    end lossConstant_D_zeta;

  partial model PartialTwoPortTransport
      "Partial element transporting fluid between two ports without storing mass or energy"
    extends Modelica_Fluid.Interfaces.PartialTwoPort(
      final port_a_exposesState=false,
      final port_b_exposesState=false);

    //Assumptions
    parameter Boolean allowFlowReversal = system.allowFlowReversal
        "allow flow reversal, false restricts to design direction (port_a -> port_b)"
      annotation(Dialog(tab="Assumptions"), Evaluate=true);

    //Initialization
    parameter Boolean compute_T = true
        "= true, if temperatures at port_a and port_b are computed" 
      annotation(Dialog(tab="Advanced"), choices(__Dymola_checkBox=true));
    parameter Medium.AbsolutePressure dp_start = 0.01*system.p_start
        "Guess value of dp = port_a.p - port_b.p" 
      annotation(Dialog(tab = "Advanced"));
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
        "Guess value of m_flow = port_a.m_flow" 
      annotation(Dialog(tab = "Advanced"));
    parameter Medium.MassFlowRate reg_m_flow_small = 0.01
        "Small mass flow rate that is used to regularize port_a_T and V_flow_a"
      annotation(Dialog(tab = "Advanced"));

    Medium.MassFlowRate m_flow(start=m_flow_start)
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
    Modelica.SIunits.VolumeFlowRate V_flow
        "Volume flow rate at inflowing port (positive when flow from port_a to port_b)";
    Modelica.SIunits.Pressure dp(start=dp_start)
        "Pressure difference between port_a and port_b (= port_a.p - port_b.p)";

    annotation (
      Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics),
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"),
      Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1})));
    Medium.Temperature port_a_T
        "Temperature close to port_a, if compute_T = true";
    Medium.Temperature port_b_T
        "Temperature close to port_b, if compute_T = true";
    Medium.ThermodynamicState port_a_state_inflow
        "Medium state close to port_a for inflowing mass flow";
    Medium.ThermodynamicState port_b_state_inflow
        "Medium state close to port_b for inflowing mass flow";
    protected
    Medium.Density port_a_d_inflow
        "Density close to port_a for inflowing mass flow";
    Medium.Density port_b_d_inflow
        "Density close to port_b for inflowing mass flow";
  equation
    // Isenthalpic state transformation (no storage and no loss of energy)
    port_a.h_outflow = inStream(port_b.h_outflow);
    port_b.h_outflow = inStream(port_a.h_outflow);

    port_a.Xi_outflow = inStream(port_b.Xi_outflow);
    port_b.Xi_outflow = inStream(port_a.Xi_outflow);

    port_a.C_outflow = inStream(port_b.C_outflow);
    port_b.C_outflow = inStream(port_a.C_outflow);

    // Medium states close to the ports when mass flows into the respective port
    port_a_state_inflow = Medium.setState_phX(port_a.p, port_b.h_outflow, port_b.Xi_outflow);
    port_b_state_inflow = Medium.setState_phX(port_b.p, port_a.h_outflow, port_a.Xi_outflow);

    // Densities close to the parts when mass flows in to the respective port
    port_a_d_inflow = Medium.density(port_a_state_inflow);
    port_b_d_inflow = Medium.density(port_b_state_inflow);

    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;

    // Design direction of mass flow rate
    m_flow = port_a.m_flow;

    // Pressure difference between ports
    dp = port_a.p - port_b.p;

    // Computation of Volume flow rate, just for plotting
    V_flow = m_flow/Modelica_Fluid.Utilities.regStep(
                     m_flow, port_a_d_inflow, port_b_d_inflow, reg_m_flow_small);

    // Computation of temperature, just for plotting
    if compute_T then
       port_a_T = Modelica_Fluid.Utilities.regStep(port_a.m_flow,
                    Medium.temperature(port_a_state_inflow),
                    Medium.temperature(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow)),
                    reg_m_flow_small);
       port_b_T = Modelica_Fluid.Utilities.regStep(port_b.m_flow,
                    Medium.temperature(port_b_state_inflow),
                    Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow)),
                    reg_m_flow_small);
    else
       port_a_T = Medium.reference_T;
       port_b_T = Medium.reference_T;
    end if;
  end PartialTwoPortTransport;
  end BaseClasses;
end Fittings;
